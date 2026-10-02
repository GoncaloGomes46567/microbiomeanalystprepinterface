from __future__ import annotations

import io
import sys
import tempfile
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from prep_core import (  # noqa: E402
    load_copy_numbers,
    load_fasta_map,
    find_group,
    process_data,
    suggest_groups_by_prefix,
    _read_table,
)

st.set_page_config(
    page_title="MicrobiomeAnalyst Prep",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    .block-container { max-width: 1180px; padding-top: 1.8rem; padding-bottom: 3rem; }
    .step { padding: .45rem .8rem; border-radius: 999px; background: #f1f3f5;
            display: inline-block; margin-right: .35rem; margin-bottom: .35rem; }
    .step.active { background: #dbeafe; font-weight: 700; }
    .muted { color: #666; font-size: .9rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

for key, value in {
    "result": None,
    "step": 1,
    "group_mode": "Prefix rules",
    "mapping": {},
    "norm_method": "Relative (TSS)",
}.items():
    st.session_state.setdefault(key, value)

st.title("🧬 MicrobiomeAnalyst Prep")
st.caption(
    "General-purpose tool to prepare abundance, taxonomy and sequence tables "
    "for the analysis workflow. No assumptions about sample names, organisms or groups."
)

steps = ["1 · Data", "2 · Groups", "3 · Review", "4 · Result"]
st.markdown(
    " ".join(
        f'<span class="step {"active" if i == st.session_state.step else ""}">{label}</span>'
        for i, label in enumerate(steps, start=1)
    ),
    unsafe_allow_html=True,
)
st.divider()


def uploader(label, help_text, types):
    return st.file_uploader(label, type=types, help=help_text)


def reset_inputs():
    for key in [
        "abundance", "taxonomy", "fasta", "copy_number_table", "result", "sample_names",
        "mapping_editor", "direct_mapping_editor",
    ]:
        st.session_state.pop(key, None)
    st.session_state.mapping = {}
    st.session_state.norm_method = "Relative (TSS)"
    st.session_state.step = 1


@st.cache_data(show_spinner=False)
def _read_table_cache(data: bytes, name: str, kind: str) -> pd.DataFrame:
    buffer = io.BytesIO(data)
    buffer.name = name  # _read_table uses the name to guess the format
    return _read_table(buffer, kind)


@st.cache_data(show_spinner=False)
def _read_fasta_cache(data: bytes):
    buffer = io.StringIO(data.decode("utf-8", errors="ignore"))
    return load_fasta_map(buffer)


@st.cache_data(show_spinner=False)
def _read_copy_numbers_cache(data: bytes, name: str):
    buffer = io.BytesIO(data)
    buffer.name = name
    return load_copy_numbers(buffer)


def read_uploaded_preview(uploaded, kind):
    return _read_table_cache(uploaded.getvalue(), uploaded.name, kind)


def read_fasta_preview(uploaded):
    return _read_fasta_cache(uploaded.getvalue())


def read_copy_numbers_preview(uploaded):
    return _read_copy_numbers_cache(uploaded.getvalue(), uploaded.name)


if st.session_state.step == 1:
    st.header("1. Upload data")
    st.write(
        "Upload the files you want to prepare. The app doesn't assume any "
        "specific set of samples or groups."
    )

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Abundance table")
        abundance = uploader(
            "Choose abundance table",
            "CSV, TSV or Excel. Should have one row per ASV and one column per sample.",
            ["tsv", "csv", "xls", "xlsx", "xlsm"],
        )
        st.caption("Required · ASV counts/abundances per sample.")

        st.subheader("Taxonomy table")
        taxonomy = uploader(
            "Choose taxonomy table",
            "CSV, TSV or Excel. The first column should identify the ASVs and another "
            "column should hold the taxonomy.",
            ["tsv", "csv", "xls", "xlsx", "xlsm"],
        )
        st.caption("Required · taxonomic classification of the ASVs.")

    with col2:
        st.subheader("ASV sequences")
        fasta = uploader(
            "Choose FASTA file",
            "Exported FASTA, or any FASTA where each sequence ID matches the ASVs.",
            ["fasta", "fa", "fna"],
        )
        st.caption("Required · sequences used to match hashes/IDs.")

        st.info(
            "💡 **Accepted formats:** CSV, TSV and Excel for tables; FASTA for sequences. "
            "The app validates the structure before processing the data."
        )

    st.subheader("Normalization")
    norm_method = st.radio(
        "Normalization method",
        ["Relative (TSS)", "16S gene copy number"],
        index=0 if st.session_state.norm_method == "Relative (TSS)" else 1,
        horizontal=True,
        help=(
            "TSS: divides each ASV by the sample's total sum (relative abundance, sums to 1). "
            "Copy number: divides each ASV by its taxon's factor in a reference table "
            "uploaded by the user — nothing is hard-coded in the app."
        ),
    )
    st.session_state.norm_method = norm_method

    copy_number_table = None
    if norm_method == "16S gene copy number":
        copy_number_table = uploader(
            "Copy number table",
            "CSV, TSV or Excel with a taxon column and a numeric factor column "
            "(e.g. 'Taxon' and 'CopyNumber'). Matching uses each ASV's most "
            "specific taxonomic level.",
            ["tsv", "csv", "xls", "xlsx", "xlsm"],
        )
        st.caption("Required for this method · no values hard-coded in the app.")

    everything_loaded = abundance and taxonomy and fasta and (
        norm_method == "Relative (TSS)" or copy_number_table
    )

    if everything_loaded:
        errors = []

        df_preview = None
        try:
            df_preview = read_uploaded_preview(abundance, "abundance")
        except Exception as exc:
            errors.append(f"Abundance table: {exc}")

        try:
            df_tax_preview = read_uploaded_preview(taxonomy, "taxonomy")
            if df_tax_preview.shape[1] == 0:
                errors.append("Taxonomy table: no taxonomy column was found.")
        except Exception as exc:
            errors.append(f"Taxonomy table: {exc}")

        fasta_collisions = []
        try:
            _, fasta_collisions = read_fasta_preview(fasta)
        except Exception as exc:
            errors.append(f"FASTA file: {exc}")

        if copy_number_table is not None:
            try:
                read_copy_numbers_preview(copy_number_table)
            except Exception as exc:
                errors.append(f"Copy number table: {exc}")

        if errors:
            for err in errors:
                st.error(err)
            st.button("Continue →", disabled=True, width="stretch")
        else:
            st.success(
                f"✓ Files selected · {df_preview.shape[0]:,} rows × "
                f"{df_preview.shape[1]:,} samples in the abundance table."
            )
            if fasta_collisions:
                st.warning(
                    f"⚠️ {len(fasta_collisions)} FASTA sequence(s) repeat a sequence already "
                    "seen under a different ID — it will be matched to the first occurrence."
                )
            with st.expander("Preview abundance table"):
                st.dataframe(df_preview.head(10), width="stretch")

            if st.button("Continue →", type="primary", width="stretch"):
                st.session_state.abundance = abundance
                st.session_state.taxonomy = taxonomy
                st.session_state.fasta = fasta
                st.session_state.copy_number_table = copy_number_table
                st.session_state.sample_names = [str(x) for x in df_preview.columns]
                st.session_state.mapping = {}
                st.session_state.step = 2
                st.rerun()
    else:
        st.button("Continue →", disabled=True, width="stretch")

elif st.session_state.step == 2:
    st.header("2. Define groups")
    st.write(
        "Choose how to assign samples to groups. You can use prefix rules "
        "or assign a group directly to each sample."
    )

    mode = st.radio(
        "Grouping method",
        ["Prefix rules", "Direct per-sample assignment"],
        index=0 if st.session_state.group_mode == "Prefix rules" else 1,
        horizontal=True,
        help=(
            "Prefix: a rule applies to every sample whose name starts with the given text. "
            "Direct: each sample is explicitly assigned a group."
        ),
    )
    st.session_state.group_mode = mode

    sample_names = st.session_state.get("sample_names", [])
    st.info(
        f"Found **{len(sample_names)} samples**. "
        "The names below come directly from the uploaded table."
    )

    if mode == "Prefix rules":
        st.markdown("### Rules")
        st.caption(
            "Generic example: prefix `A` could map to group `Group 1`. "
            "There are no predefined groups. If prefixes overlap, the most "
            "specific rule is used first."
        )

        col_suggest, _ = st.columns([1, 2])
        with col_suggest:
            if st.button("💡 Suggest from sample names", width="stretch"):
                suggestion = suggest_groups_by_prefix(sample_names)
                if suggestion:
                    st.session_state.mapping = suggestion
                    st.session_state.pop("mapping_editor", None)
                    st.rerun()
                else:
                    st.info(
                        "No repeated prefixes were found among the samples "
                        "(every name looks unique). Define the groups manually."
                    )

        current = st.session_state.mapping
        mapping_df = pd.DataFrame(
            [{"Prefix": k, "Group": v} for k, v in current.items()],
            columns=["Prefix", "Group"],
        )
        if mapping_df.empty:
            mapping_df = pd.DataFrame([{"Prefix": "", "Group": ""}])

        edited = st.data_editor(
            mapping_df,
            num_rows="dynamic",
            width="stretch",
            hide_index=True,
            column_config={
                "Prefix": st.column_config.TextColumn(
                    "Sample prefix", required=False,
                    help="Leading text of the sample name."
                ),
                "Group": st.column_config.TextColumn(
                    "Group name", required=False,
                    help="Free-form name chosen by the user."
                ),
            },
            key="mapping_editor",
        )

        if sample_names:
            rules = []
            for _, row in edited.iterrows():
                p = str(row.get("Prefix", "")).strip()
                g = str(row.get("Group", "")).strip()
                if p and g and p.lower() != "nan" and g.lower() != "nan":
                    rules.append((p, g))

            if rules:
                def preview_group(name):
                    for p, g in sorted(rules, key=lambda x: len(x[0]), reverse=True):
                        if name.startswith(p):
                            return g
                    return "No group"

                preview = pd.DataFrame({
                    "Sample": sample_names,
                    "Detected group": [preview_group(n) for n in sample_names],
                })
                with st.expander("See how the rules will be applied"):
                    st.dataframe(preview, width="stretch", hide_index=True)
                    missing = int((preview["Detected group"] == "No group").sum())
                    if missing:
                        st.warning(f"{missing} sample(s) still don't match any rule.")

    else:
        st.markdown("### Sample assignment")

        col_suggest, _ = st.columns([1, 2])
        with col_suggest:
            if st.button("💡 Suggest from sample names", width="stretch"):
                suggested_prefixes = suggest_groups_by_prefix(sample_names)
                if suggested_prefixes:
                    st.session_state.mapping = {
                        name: group
                        for name in sample_names
                        if (group := find_group(name, suggested_prefixes)) != "Unknown"
                    }
                    st.session_state.pop("direct_mapping_editor", None)
                    st.rerun()
                else:
                    st.info(
                        "No repeated prefixes were found among the samples "
                        "(every name looks unique). Assign the groups manually."
                    )

        direct_df = pd.DataFrame({
            "Sample": sample_names,
            "Group": [
                st.session_state.mapping.get(name, "") for name in sample_names
            ],
        })
        edited = st.data_editor(
            direct_df,
            width="stretch",
            hide_index=True,
            disabled=["Sample"],
            column_config={
                "Sample": st.column_config.TextColumn("Sample"),
                "Group": st.column_config.TextColumn(
                    "Group", help="Enter the group for each sample."
                ),
            },
            key="direct_mapping_editor",
        )
        st.caption("You can leave a sample without a group; it will be marked as `Unknown`.")

    c1, c2 = st.columns(2)
    with c1:
        if st.button("← Back", width="stretch"):
            st.session_state.step = 1
            st.rerun()
    with c2:
        if st.button("Continue →", type="primary", width="stretch"):
            clean = {}
            if mode == "Prefix rules":
                for _, row in edited.iterrows():
                    prefix = str(row.get("Prefix", "")).strip()
                    group = str(row.get("Group", "")).strip()
                    if prefix and group and prefix.lower() != "nan" and group.lower() != "nan":
                        clean[prefix] = group
            else:
                for _, row in edited.iterrows():
                    sample = str(row.get("Sample", "")).strip()
                    group = str(row.get("Group", "")).strip()
                    if sample and group and group.lower() != "nan":
                        clean[sample] = group

            st.session_state.mapping = clean
            st.session_state.step = 3
            st.rerun()

elif st.session_state.step == 3:
    st.header("3. Review and run")
    st.write("Confirm the data and rules before starting the processing.")

    st.subheader("Files")
    for label, key in [
        ("Abundance table", "abundance"),
        ("Taxonomy table", "taxonomy"),
        ("FASTA sequences", "fasta"),
    ]:
        st.write(f"**{label}:** `{st.session_state[key].name}`")
    if st.session_state.get("copy_number_table") is not None:
        st.write(f"**Copy number table:** `{st.session_state.copy_number_table.name}`")

    st.subheader("Normalization")
    st.write(f"**Method:** {st.session_state.norm_method}")

    st.subheader("Group configuration")
    if st.session_state.mapping:
        st.dataframe(
            pd.DataFrame([
                {"Rule / sample": k, "Group": v}
                for k, v in st.session_state.mapping.items()
            ]),
            width="stretch",
            hide_index=True,
        )
    else:
        st.warning(
            "No group rule was defined. Samples will be processed and will "
            "appear as `Unknown` in the metadata."
        )

    norm_label = "TSS" if st.session_state.norm_method == "Relative (TSS)" else "copy number"
    st.info(
        f"Processing will generate raw abundance, normalized abundance "
        f"({norm_label}), taxonomy and metadata files."
    )

    c1, c2 = st.columns(2)
    with c1:
        if st.button("← Change groups", width="stretch"):
            st.session_state.step = 2
            st.rerun()
    with c2:
        if st.button("▶ Run processing", type="primary", width="stretch"):
            with st.spinner("Preparing data…"):
                try:
                    with tempfile.TemporaryDirectory() as tmp:
                        copy_number_table = st.session_state.get("copy_number_table")
                        for uploaded in (
                            st.session_state.abundance,
                            st.session_state.taxonomy,
                            st.session_state.fasta,
                            copy_number_table,
                        ):
                            if uploaded is not None:
                                uploaded.seek(0)
                        method = "tss" if st.session_state.norm_method == "Relative (TSS)" else "copy_number"
                        result = process_data(
                            st.session_state.abundance,
                            st.session_state.fasta,
                            st.session_state.taxonomy,
                            st.session_state.mapping,
                            tmp,
                            normalization_method=method,
                            copy_number_table=copy_number_table,
                        )
                        result["file_bytes"] = {
                            path.name: Path(path).read_bytes() for path in result["files"]
                        }
                        st.session_state.result = result
                    st.session_state.step = 4
                    st.rerun()
                except Exception as exc:
                    st.error(f"Processing could not be completed:\n\n{exc}")

elif st.session_state.step == 4:
    result = st.session_state.result
    st.header("4. Result")

    if result is None:
        st.warning("No result yet.")
        if st.button("Back to start"):
            reset_inputs()
            st.rerun()
        st.stop()

    if result["unknown_samples"]:
        st.warning(
            f"{len(result['unknown_samples'])} sample(s) ended up without a group. "
            "They will be marked as `Unknown` in metadata.txt."
        )
        with st.expander("See samples without a group"):
            st.write(", ".join(result["unknown_samples"]))
    else:
        st.success("✓ All samples were assigned to a group.")

    if result["fasta_unmatched"]:
        with st.expander(
            f"ℹ️ {len(result['fasta_unmatched'])} ASV(s) kept their original ID "
            "(no match in the FASTA)"
        ):
            st.caption(
                "This is expected if the table already uses readable IDs instead of "
                "sequence hashes. If you weren't expecting this, check that you "
                "uploaded the right FASTA file."
            )
            st.write(", ".join(result["fasta_unmatched"][:100]))

    if result["fasta_collisions"]:
        st.warning(
            f"⚠️ {len(result['fasta_collisions'])} FASTA sequence(s) were identical "
            "to a sequence already processed and were skipped."
        )

    if result["copy_number_unmatched"]:
        with st.expander(
            f"ℹ️ {len(result['copy_number_unmatched'])} ASV(s) without a match "
            "in the copy number table (factor 1.0 applied)"
        ):
            st.caption(
                "The most specific taxon for these ASVs wasn't found in the uploaded "
                "reference table. For a more accurate factor, add that taxon to the "
                "table and run the processing again."
            )
            st.write(", ".join(result["copy_number_unmatched"][:100]))

    m1, m2, m3 = st.columns(3)
    m1.metric("ASVs", result["asvs"])
    m2.metric("Samples", result["samples"])
    m3.metric("Groups", result["groups"])

    st.subheader("Generated files")
    norm_description = (
        "TSS-normalized table (relative abundance)"
        if result["normalization_method"] == "tss"
        else "Table normalized by 16S gene copy number"
    )
    files_df = pd.DataFrame({
        "File": list(result["file_bytes"].keys()),
        "Description": [
            "Raw abundance table",
            norm_description,
            "Formatted taxonomy",
            "Metadata and groups",
        ],
    })
    st.dataframe(files_df, width="stretch", hide_index=True)

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for filename, content in result["file_bytes"].items():
            zf.writestr(filename, content)
    zip_buffer.seek(0)

    st.download_button(
        "⬇️ Download all results (ZIP)",
        data=zip_buffer.getvalue(),
        file_name="phase2_results.zip",
        mime="application/zip",
        type="primary",
        width="stretch",
    )

    st.subheader("Preview")
    tab1, tab2, tab3 = st.tabs(["Metadata", "Taxonomy", "Abundances"])
    with tab1:
        st.dataframe(result["metadata"].head(50), width="stretch")
    with tab2:
        st.dataframe(result["taxonomy"].head(50), width="stretch")
    with tab3:
        st.dataframe(result["otu_table"].iloc[:50, :50], width="stretch")

    if st.button("↻ Prepare new data", width="stretch"):
        reset_inputs()
        st.rerun()
