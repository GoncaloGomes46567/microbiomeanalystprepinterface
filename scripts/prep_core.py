from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Tuple, Union

import pandas as pd

PathLike = Union[str, Path]

NIVEIS_TAXONOMICOS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PREFIXO_PARA_NIVEL = {
    "k__": "Kingdom", "d__": "Kingdom", "p__": "Phylum", "c__": "Class",
    "o__": "Order", "f__": "Family", "g__": "Genus", "s__": "Species",
}


def _open_text(path_or_file):
    """Return a text stream, either for a path or an already opened upload.

    We don't close UploadedFile objects here because Streamlit might need
    to reuse them within the same interaction.
    """
    if hasattr(path_or_file, "read"):
        return path_or_file, False
    return open(path_or_file, "r", encoding="utf-8"), True


def carregar_mapa_fasta(caminho_fasta) -> Tuple[Dict[str, str], List[str]]:
    """Build a mapping from md5(sequence) -> ASV ID based on a FASTA file.

    Also returns a list of IDs whose sequence collides with that of a previous ID
    (same sequence, different header), so the user can be warned instead of the collision going unnoticed.
    """
    mapa: Dict[str, str] = {}
    colisoes: List[str] = []

    def registar(id_atual, partes):
        if not id_atual or not partes:
            return
        chave = hashlib.md5("".join(partes).encode()).hexdigest()
        if chave in mapa and mapa[chave] != id_atual:
            colisoes.append(id_atual)
        else:
            mapa[chave] = id_atual

    f, should_close = _open_text(caminho_fasta)
    try:
        id_atual, partes = None, []
        for linha in f:
            linha = linha.decode("utf-8") if isinstance(linha, bytes) else linha
            linha = linha.strip()
            if not linha:
                continue
            if linha.startswith(">"):
                registar(id_atual, partes)
                id_atual = linha[1:].split()[0]
                partes = []
            else:
                partes.append(linha)
        registar(id_atual, partes)
    finally:
        if should_close:
            f.close()

    if not mapa:
        raise ValueError("The FASTA file does not contain valid sequences.")
    return mapa, colisoes


def limpar_taxonomia(texto_tax):
    """Separate a string of the type 'k__Bacteria;p__Firmicutes;...' by its levels."""
    resultado = {nivel: "" for nivel in NIVEIS_TAXONOMICOS}
    if pd.isna(texto_tax):
        return resultado
    for parte in str(texto_tax).split(";"):
        parte = parte.strip()
        for prefixo, nivel in PREFIXO_PARA_NIVEL.items():
            if parte.startswith(prefixo):
                resultado[nivel] = parte[len(prefixo):].strip()
    return resultado


def formatar_taxonomia(df_tax: pd.DataFrame) -> pd.DataFrame:
    """Return a table Kingdom..Species based on two possible formats:

    - uma única coluna com strings 'k__x;p__y;...' (export QIIME2/greengenes); ou
    - colunas já separadas com os nomes dos níveis (ex.: export do MicrobiomeAnalyst).

   if neither format is detected, the function will return a DataFrame with empty strings for all levels.
    """
    colunas_por_nome = {str(c).strip().lower(): c for c in df_tax.columns}
    tem_colunas_separadas = any(nivel.lower() in colunas_por_nome for nivel in NIVEIS_TAXONOMICOS)

    if tem_colunas_separadas:
        tax = pd.DataFrame(index=df_tax.index)
        for nivel in NIVEIS_TAXONOMICOS:
            coluna = colunas_por_nome.get(nivel.lower())
            tax[nivel] = df_tax[coluna].fillna("") if coluna is not None else ""
        return tax

    return df_tax[df_tax.columns[0]].apply(limpar_taxonomia).apply(pd.Series)


def descobrir_grupo(nome_amostra: str, dicionario_grupos: Mapping[str, str]) -> str:
    """Determine the group of a sample based on its name and a mapping of prefixes to groups."""
    for prefixo in sorted(dicionario_grupos.keys(), key=len, reverse=True):
        if nome_amostra.startswith(prefixo):
            return dicionario_grupos[prefixo]
    return "Unknown"


def sugerir_grupos_por_prefixo(nomes_amostra: Iterable[str]) -> Dict[str, str]:
    """Propose rules → group only based on the names of the loaded samples.

    There is no predefined group (neither "Frutalose", nor "Inoculum", nor
    anything specific to a particular study): the idea is to look at the part
    of each name that normally identifies the replica (a number, or a
    number followed by a dot/hyphen/underscore at the end) and remove that
    parte. O que sobra é usado como prefixo do grupo.

    Examples:
        "A.CN.1", "A.CN.2"      -> prefixo "A.CN" (2 amostras, sugerido)
        "V1", "V2", "V3"        -> prefixo "V" (3 amostras, sugerido)
        "amostra_unica"         -> não é sugerido (não repete com mais nenhuma)

    Only suggest a prefix when at least two samples share it —
    otherwise there are no replicas to group and the suggestion would
    not be helpful. The name of the suggested group is the prefix itself; the
    user can always edit it in the table before continuing.
    """
    from collections import Counter

    prefixo_por_amostra = {}
    for nome in nomes_amostra:
        nome = str(nome)
        nucleo = re.sub(r"[\s._-]*\d+$", "", nome).strip()
        prefixo_por_amostra[nome] = nucleo or nome

    contagem = Counter(prefixo_por_amostra.values())
    prefixos_com_replicas = [p for p, n in contagem.items() if n >= 2 and p]

 
    return {prefixo: prefixo for prefixo in sorted(prefixos_com_replicas, key=len, reverse=True)}


def _detetar_separador(primeira_linha: str) -> str:
    """Try to detect the separator used in the first line."""
    if primeira_linha.count(",") > primeira_linha.count("\t"):
        return ","
    return "\t"


def _ler_primeiras_linhas(file_or_path, n: int = 2) -> List[str]:
    """Read `n` lines from a file or upload, without consuming the stream."""
    if hasattr(file_or_path, "read"):
        file_or_path.seek(0)
        conteudo = file_or_path.read()
        file_or_path.seek(0)
    else:
        with open(file_or_path, "rb") as fh:
            conteudo = fh.read()

    if isinstance(conteudo, bytes):
        conteudo = conteudo.decode("utf-8", errors="ignore")
    return conteudo.splitlines()[:n]


def _read_table(file_or_path, kind: str) -> pd.DataFrame:
    """Read CSV/TSV/XLS/XLSX, detecting the separator by its content instead of assuming it by the extension.

    Also recognizes the comment header that QIIME2 places before the
    column header line in exported abundance tables.
    """
    nome = getattr(file_or_path, "name", str(file_or_path)).lower()
    sufixo = Path(nome).suffix.lower()

    if sufixo in {".xls", ".xlsx", ".xlsm"}:
        return pd.read_excel(file_or_path, index_col=0)

    linhas = _ler_primeiras_linhas(file_or_path)
    primeira_linha = linhas[0] if linhas else ""

    separador = "," if sufixo == ".csv" else _detetar_separador(primeira_linha)

   
    ignora_primeira = (
        kind == "abundance"
        and primeira_linha.startswith("#")
        and not primeira_linha.lower().lstrip("#").startswith(("otu id", "name"))
    )

    if hasattr(file_or_path, "seek"):
        file_or_path.seek(0)

    return pd.read_csv(
        file_or_path,
        sep=separador,
        skiprows=1 if ignora_primeira else 0,
        index_col=0,
    )


def validate_inputs(df_otu: pd.DataFrame, df_tax: pd.DataFrame, fasta_map: Mapping[str, str]):
    """Validate the input dataframes."""
    errors = []
    if df_otu.empty:
        errors.append("Abundance table is empty.")
    if df_tax.empty:
        errors.append("Taxonomy table is empty.")
    if df_otu.index.duplicated().any():
        errors.append("Abundance table contains duplicate ASV IDs.")
    if df_tax.index.duplicated().any():
        errors.append("Taxonomy table contains duplicate ASV IDs.")
    if df_otu.shape[1] == 0:
        errors.append("No samples found in the abundance table.")

    try:
        numeric = df_otu.apply(pd.to_numeric)
        if numeric.isna().any().any():
            errors.append("The abundance table contains non-numeric values.")
    except Exception:
        errors.append("Could not interpret the abundance table as numeric values.")

    if not fasta_map:
        errors.append("No sequences found in the FASTA file.")
    return errors


def process_data(
    tabela_input,
    fasta,
    taxonomia,
    mapeamento: Mapping[str, str],
    pasta_saida: PathLike,
) -> dict:
    """Run the complete pipeline for Phase 2 and return a summary of the results."""
    output_dir = Path(pasta_saida)
    output_dir.mkdir(parents=True, exist_ok=True)

    df_otu = _read_table(tabela_input, "abundance")
    if hasattr(tabela_input, "seek"):
        tabela_input.seek(0)

    df_tax = _read_table(taxonomia, "taxonomy")
    if hasattr(taxonomia, "seek"):
        taxonomia.seek(0)

    mapa_ids, sequencias_duplicadas = carregar_mapa_fasta(fasta)
    errors = validate_inputs(df_otu, df_tax, mapa_ids)
    if errors:
        raise ValueError("\n".join(errors))

    ids_sem_correspondencia = [str(idx) for idx in df_otu.index if str(idx) not in mapa_ids]

    df_otu = df_otu.apply(pd.to_numeric)
    df_otu.index = [mapa_ids.get(str(idx), str(idx)) for idx in df_otu.index]
    df_otu = df_otu.astype(int)
    df_otu.index.name = "#NAME"
    otu_path = output_dir / "otu_table.txt"
    df_otu.to_csv(otu_path, sep="\t")

    col_sums = df_otu.sum(axis=0)
    if (col_sums <= 0).any():
        bad = ", ".join(map(str, col_sums[col_sums <= 0].index))
        raise ValueError(f"Groups with zero counts: {bad}")

    otu_tss = df_otu.div(col_sums, axis=1)
    desvio_max = float((otu_tss.sum(axis=0) - 1.0).abs().max())
    if desvio_max >= 1e-9:
        raise ValueError(f"TSS validation failed (maximum deviation: {desvio_max:.2e}).")
    otu_tss.index.name = "#NAME"
    tss_path = output_dir / "otu_table_tss.txt"
    otu_tss.to_csv(tss_path, sep="\t", float_format="%.15g")

    if df_tax.shape[1] == 0:
        raise ValueError("The taxonomy table does not contain any taxonomy columns.")
    df_tax.index = [mapa_ids.get(str(idx), str(idx)) for idx in df_tax.index]
    tax_limpa = formatar_taxonomia(df_tax)
    tax_limpa = tax_limpa.reindex(df_otu.index).fillna("")
    tax_limpa.index.name = "#TAXONOMY"
    tax_path = output_dir / "taxonomy.txt"
    tax_limpa.to_csv(tax_path, sep="\t")

    lista_grupos = [descobrir_grupo(str(a), mapeamento) for a in df_otu.columns]
    desconhecidos = [s for s, g in zip(df_otu.columns, lista_grupos) if g == "Unknown"]
    df_meta = pd.DataFrame({"Group": lista_grupos}, index=df_otu.columns)
    df_meta.index.name = "#NAME"
    meta_path = output_dir / "metadata.txt"
    df_meta.to_csv(meta_path, sep="\t")

    return {
        "files": [otu_path, tss_path, tax_path, meta_path],
        "asvs": int(df_otu.shape[0]),
        "samples": int(df_otu.shape[1]),
        "groups": int(df_meta["Group"].nunique()),
        "unknown_samples": desconhecidos,
        "tss_max_deviation": desvio_max,
        "fasta_sem_correspondencia": ids_sem_correspondencia,
        "fasta_sequencias_duplicadas": sequencias_duplicadas,
        "otu_table": df_otu,
        "taxonomy": tax_limpa,
        "metadata": df_meta,
    }


def load_mapping(path: PathLike) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return {str(k): str(v) for k, v in data.items()}
