from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Tuple, Union

import pandas as pd

PathLike = Union[str, Path]

# SILVA/GTDB use d__ instead of k__ for kingdom
TAXONOMY_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PREFIX_TO_LEVEL = {
    "k__": "Kingdom", "d__": "Kingdom", "p__": "Phylum", "c__": "Class",
    "o__": "Order", "f__": "Family", "g__": "Genus", "s__": "Species",
}


def _open_text(path_or_file):
    if hasattr(path_or_file, "read"):
        return path_or_file, False
    return open(path_or_file, "r", encoding="utf-8"), True


def load_fasta_map(fasta_path) -> Tuple[Dict[str, str], List[str]]:
    seq_map: Dict[str, str] = {}
    collisions: List[str] = []

    def register(current_id, parts):
        if not current_id or not parts:
            return
        key = hashlib.md5("".join(parts).encode()).hexdigest()
        if key in seq_map and seq_map[key] != current_id:
            collisions.append(current_id)
        else:
            seq_map[key] = current_id

    f, should_close = _open_text(fasta_path)
    try:
        current_id, parts = None, []
        for line in f:
            line = line.decode("utf-8") if isinstance(line, bytes) else line
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                register(current_id, parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        register(current_id, parts)
    finally:
        if should_close:
            f.close()

    if not seq_map:
        raise ValueError("FASTA file has no valid sequences.")
    return seq_map, collisions


def clean_taxonomy(raw_taxon):
    result = {level: "" for level in TAXONOMY_LEVELS}
    if pd.isna(raw_taxon):
        return result
    for part in str(raw_taxon).split(";"):
        part = part.strip()
        for prefix, level in PREFIX_TO_LEVEL.items():
            if part.startswith(prefix):
                result[level] = part[len(prefix):].strip()
    return result


def format_taxonomy(tax_df: pd.DataFrame) -> pd.DataFrame:
    columns_by_name = {str(c).strip().lower(): c for c in tax_df.columns}
    has_separate_columns = any(level.lower() in columns_by_name for level in TAXONOMY_LEVELS)

    if has_separate_columns:
        tax = pd.DataFrame(index=tax_df.index)
        for level in TAXONOMY_LEVELS:
            col = columns_by_name.get(level.lower())
            tax[level] = tax_df[col].fillna("") if col is not None else ""
        return tax

    return tax_df[tax_df.columns[0]].apply(clean_taxonomy).apply(pd.Series)


def _last_assigned_level(tax_row: pd.Series) -> str:
    for level in reversed(TAXONOMY_LEVELS):
        value = str(tax_row.get(level, "")).strip()
        if value and value.lower() != "nan":
            return value
    return ""


def load_copy_numbers(file_or_path) -> Dict[str, float]:
    df = _read_table(file_or_path, "copynumber")
    df = df.reset_index()  # the taxon column may have been read as the index
    columns = {str(c).strip().lower(): c for c in df.columns}

    def find(patterns):
        for name_lower, name in columns.items():
            if any(p in name_lower for p in patterns):
                return name
        return None

    taxon_col = find(["taxon", "taxa", "name"])
    copy_col = find(["copy", "cn"])
    if taxon_col is None or copy_col is None:
        raise ValueError(
            "Could not find taxon/copy-number columns in the reference table "
            "(expected something like 'Taxon' and 'CopyNumber')."
        )

    copy_map: Dict[str, float] = {}
    for taxon, copy_number in zip(df[taxon_col], df[copy_col]):
        taxon = str(taxon).strip()
        if not taxon or taxon.lower() == "nan":
            continue
        try:
            copy_map[taxon.lower()] = float(copy_number)
        except (TypeError, ValueError):
            continue

    if not copy_map:
        raise ValueError("Copy-number table is empty or has no valid numeric values.")
    return copy_map


def normalize_by_copy_number(
    otu_df: pd.DataFrame,
    tax_df: pd.DataFrame,
    copy_map: Mapping[str, float],
) -> Tuple[pd.DataFrame, List[str]]:
    unmatched: List[str] = []
    factors = []

    for asv in otu_df.index:
        taxon = _last_assigned_level(tax_df.loc[asv]) if asv in tax_df.index else ""
        taxon_lower = taxon.lower()
        factor = copy_map.get(taxon_lower)

        if factor is None and taxon_lower:
            for name, value in copy_map.items():
                if taxon_lower in name or name in taxon_lower:
                    factor = value
                    break

        if factor is None:
            factor = 1.0
            unmatched.append(str(asv))

        factors.append(factor)

    factors_series = pd.Series(factors, index=otu_df.index)
    normalized = otu_df.div(factors_series, axis=0).clip(lower=0)
    return normalized, unmatched


def find_group(sample_name: str, group_rules: Mapping[str, str]) -> str:
    for prefix in sorted(group_rules.keys(), key=len, reverse=True):
        if sample_name.startswith(prefix):
            return group_rules[prefix]
    return "Unknown"


def suggest_groups_by_prefix(sample_names: Iterable[str]) -> Dict[str, str]:
    from collections import Counter

    prefix_by_sample = {}
    for name in sample_names:
        name = str(name)
        core = re.sub(r"[\s._-]*\d+$", "", name).strip()
        prefix_by_sample[name] = core or name

    counts = Counter(prefix_by_sample.values())
    repeated = [p for p, n in counts.items() if n >= 2 and p]

    return {prefix: prefix for prefix in sorted(repeated, key=len, reverse=True)}


def _detect_separator(first_line: str) -> str:
    return "," if first_line.count(",") > first_line.count("\t") else "\t"


def _read_first_lines(file_or_path, n: int = 2) -> List[str]:
    if hasattr(file_or_path, "read"):
        file_or_path.seek(0)
        content = file_or_path.read()
        file_or_path.seek(0)
    else:
        with open(file_or_path, "rb") as fh:
            content = fh.read()

    if isinstance(content, bytes):
        content = content.decode("utf-8", errors="ignore")
    return content.splitlines()[:n]


def _read_table(file_or_path, kind: str) -> pd.DataFrame:
    name = getattr(file_or_path, "name", str(file_or_path)).lower()
    suffix = Path(name).suffix.lower()

    if suffix in {".xls", ".xlsx", ".xlsm"}:
        return pd.read_excel(file_or_path, index_col=0)

    lines = _read_first_lines(file_or_path)
    first_line = lines[0] if lines else ""

    sep = "," if suffix == ".csv" else _detect_separator(first_line)

    # QIIME2 exports have a comment line before the real header
    skip_first = (
        kind == "abundance"
        and first_line.startswith("#")
        and not first_line.lower().lstrip("#").startswith(("otu id", "name"))
    )

    if hasattr(file_or_path, "seek"):
        file_or_path.seek(0)

    return pd.read_csv(
        file_or_path,
        sep=sep,
        skiprows=1 if skip_first else 0,
        index_col=0,
    )


def validate_inputs(otu_df: pd.DataFrame, tax_df: pd.DataFrame, fasta_map: Mapping[str, str]):
    errors = []
    if otu_df.empty:
        errors.append("Abundance table is empty.")
    if tax_df.empty:
        errors.append("Taxonomy table is empty.")
    if otu_df.index.duplicated().any():
        errors.append("Abundance table has duplicate ASV IDs.")
    if tax_df.index.duplicated().any():
        errors.append("Taxonomy table has duplicate ASV IDs.")
    if otu_df.shape[1] == 0:
        errors.append("No samples found in the abundance table.")

    try:
        numeric = otu_df.apply(pd.to_numeric)
        if numeric.isna().any().any():
            errors.append("Abundance table contains non-numeric values.")
    except Exception:
        errors.append("Could not parse the abundance table as numeric values.")

    if not fasta_map:
        errors.append("No sequences found in the FASTA file.")
    return errors


def process_data(
    input_table,
    fasta,
    taxonomy,
    group_mapping: Mapping[str, str],
    output_folder: PathLike,
    normalization_method: str = "tss",
    copy_number_table=None,
) -> dict:
    if normalization_method not in {"tss", "copy_number"}:
        raise ValueError(f"Unknown normalization method: {normalization_method}")
    if normalization_method == "copy_number" and copy_number_table is None:
        raise ValueError("A copy-number reference table is required for this normalization method.")

    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    otu_df = _read_table(input_table, "abundance")
    if hasattr(input_table, "seek"):
        input_table.seek(0)

    tax_df = _read_table(taxonomy, "taxonomy")
    if hasattr(taxonomy, "seek"):
        taxonomy.seek(0)

    seq_map, fasta_collisions = load_fasta_map(fasta)
    errors = validate_inputs(otu_df, tax_df, seq_map)
    if errors:
        raise ValueError("\n".join(errors))

    fasta_unmatched = [str(idx) for idx in otu_df.index if str(idx) not in seq_map]

    otu_df = otu_df.apply(pd.to_numeric)
    otu_df.index = [seq_map.get(str(idx), str(idx)) for idx in otu_df.index]
    otu_df = otu_df.astype(int)
    otu_df.index.name = "#NAME"
    otu_path = output_dir / "otu_table.txt"
    otu_df.to_csv(otu_path, sep="\t")

    col_sums = otu_df.sum(axis=0)
    if (col_sums <= 0).any():
        bad = ", ".join(map(str, col_sums[col_sums <= 0].index))
        raise ValueError(f"Samples with zero total counts: {bad}")

    if tax_df.shape[1] == 0:
        raise ValueError("Taxonomy table has no taxonomy column.")
    tax_df.index = [seq_map.get(str(idx), str(idx)) for idx in tax_df.index]
    clean_tax = format_taxonomy(tax_df)
    clean_tax = clean_tax.reindex(otu_df.index).fillna("")
    clean_tax.index.name = "#TAXONOMY"
    tax_path = output_dir / "taxonomy.txt"
    clean_tax.to_csv(tax_path, sep="\t")

    copy_number_unmatched: List[str] = []
    tss_max_deviation = None
    if normalization_method == "tss":
        norm_df = otu_df.div(col_sums, axis=1)
        tss_max_deviation = float((norm_df.sum(axis=0) - 1.0).abs().max())
        if tss_max_deviation >= 1e-9:
            raise ValueError(f"TSS validation failed (max deviation: {tss_max_deviation:.2e}).")
        norm_path = output_dir / "otu_table_tss.txt"
    else:
        copy_map = load_copy_numbers(copy_number_table)
        norm_df, copy_number_unmatched = normalize_by_copy_number(otu_df, clean_tax, copy_map)
        norm_path = output_dir / "otu_table_copynumber.txt"

    norm_df.index.name = "#NAME"
    norm_df.to_csv(norm_path, sep="\t", float_format="%.15g")

    groups = [find_group(str(s), group_mapping) for s in otu_df.columns]
    unknown_samples = [s for s, g in zip(otu_df.columns, groups) if g == "Unknown"]
    metadata_df = pd.DataFrame({"Group": groups}, index=otu_df.columns)
    metadata_df.index.name = "#NAME"
    metadata_path = output_dir / "metadata.txt"
    metadata_df.to_csv(metadata_path, sep="\t")

    return {
        "files": [otu_path, norm_path, tax_path, metadata_path],
        "asvs": int(otu_df.shape[0]),
        "samples": int(otu_df.shape[1]),
        "groups": int(metadata_df["Group"].nunique()),
        "unknown_samples": unknown_samples,
        "normalization_method": normalization_method,
        "tss_max_deviation": tss_max_deviation,
        "copy_number_unmatched": copy_number_unmatched,
        "fasta_unmatched": fasta_unmatched,
        "fasta_collisions": fasta_collisions,
        "otu_table": otu_df,
        "normalized_otu_table": norm_df,
        "taxonomy": clean_tax,
        "metadata": metadata_df,
    }


def load_mapping(path: PathLike) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return {str(k): str(v) for k, v in data.items()}
