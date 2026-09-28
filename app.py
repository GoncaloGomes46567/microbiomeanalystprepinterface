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
    carregar_mapa_fasta,
    descobrir_grupo,
    process_data,
    sugerir_grupos_por_prefixo,
    _read_table,
)

st.set_page_config(
    page_title="MicrobiomeAnalyst Data Prep",
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
    "group_mode": "Rules by prefix",
    "mapping": {},
}.items():
    st.session_state.setdefault(key, value)

st.title("🧬 MicrobiomeAnalyst Data Prep")
st.caption(
    "Tool to prepare files to upload to MicrobiomeAnalyst. "
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
        "abundance", "taxonomy", "fasta", "result", "sample_names",
        "mapping_editor", "direct_mapping_editor",
    ]:
        st.session_state.pop(key, None)
    st.session_state.mapping = {}
    st.session_state.step = 1


@st.cache_data(show_spinner=False)
def _ler_tabela_cache(dados: bytes, nome: str, kind: str) -> pd.DataFrame:
    buffer = io.BytesIO(dados)
    buffer.name = nome  
    return _read_table(buffer, kind)


@st.cache_data(show_spinner=False)
def _ler_fasta_cache(dados: bytes):
    buffer = io.StringIO(dados.decode("utf-8", errors="ignore"))
    return carregar_mapa_fasta(buffer)


def read_uploaded_preview(uploaded, kind):
    return _ler_tabela_cache(uploaded.getvalue(), uploaded.name, kind)


def read_fasta_preview(uploaded):
    return _ler_fasta_cache(uploaded.getvalue())


if st.session_state.step == 1:
    st.header("1. Data upload")
    st.write(
        "Upload the files that you want to process. The app will validate the structure before processing the data."
    )

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Abundance table")
        abundance = uploader(
            "Choose abundance table",
            "CSV, TSV or Excel. Should contain one row per ASV and one column per sample.",
            ["tsv", "csv", "xls", "xlsx", "xlsm"],
        )
        st.caption("Obrigatório · contagens/abundâncias das ASVs por amostra.")

        st.subheader("Tabela de taxonomia")
        taxonomy = uploader(
            "Choose taxonomy table",
            "CSV, TSV or Excel. The first column should identify the ASVs and another column should contain the taxonomy.",
            ["tsv", "csv", "xls", "xlsx", "xlsm"],
        )
        st.caption("Obrigatório · identificação/classificação taxonómica das ASVs.")

    with col2:
        st.subheader("Sequências ASV")
        fasta = uploader(
            "Choose FASTA file",
            "Exported FASTA or another FASTA where each sequence identifier corresponds to the ASVs.",
            ["fasta", "fa", "fna"],
        )
        st.caption("Required · sequences used to associate hashes/IDs.")

        st.info(
            "💡 **Accepted formats:** CSV, TSV and Excel for tables; FASTA for sequences. "
            "The app validates the structure before processing the data."
        )

    if abundance and taxonomy and fasta:
        erros = []

        df_preview = None
        try:
            df_preview = read_uploaded_preview(abundance, "abundance")
        except Exception as exc:
            erros.append(f"Abundance table: {exc}")

        try:
            df_tax_preview = read_uploaded_preview(taxonomy, "taxonomy")
            if df_tax_preview.shape[1] == 0:
                erros.append("Taxonomy table: no column with taxonomic information found.")
        except Exception as exc:
            erros.append(f"Taxonomy table: {exc}")

        colisoes_fasta = []
        try:
            _, colisoes_fasta = read_fasta_preview(fasta)
        except Exception as exc:
            erros.append(f"FASTA file: {exc}")

        if erros:
            for erro in erros:
                st.error(erro)
            st.button("Continuar →", disabled=True, width="stretch")
        else:
            st.success(
                f"✓ Files selected · {df_preview.shape[0]:,} rows × "
                f"{df_preview.shape[1]:,} samples in the abundance table."
            )
            if colisoes_fasta:
                st.warning(
                    f"!!! {len(colisoes_fasta)} sequence(s) in the FASTA file repeat a sequence already seen "
                    "with a different identifier — it will be associated with the first occurrence."
                )
            with st.expander("Pré-visualizar tabela de abundâncias"):
                st.dataframe(df_preview.head(10), width="stretch")

            if st.button("Continuar →", type="primary", width="stretch"):
                st.session_state.abundance = abundance
                st.session_state.taxonomy = taxonomy
                st.session_state.fasta = fasta
                st.session_state.sample_names = [str(x) for x in df_preview.columns]
                st.session_state.mapping = {}
                st.session_state.step = 2
                st.rerun()
    else:
        st.button("Continuar →", disabled=True, width="stretch")

elif st.session_state.step == 2:
    st.header("2. Definir grupos")
    st.write(
        "Choose how you want to associate the samples with the groups. "
        "You can use rules by prefix or assign a group directly to each sample."
    )

    mode = st.radio(
        "Grouping method",
        ["Rules by prefix", "Direct assignment by sample"],
        index=0 if st.session_state.group_mode == "Regras por prefixo" else 1,
        horizontal=True,
        help=(
            "Prefix: a rule applies to all samples whose name starts with the indicated text. "
            "Direct: each sample is explicitly assigned a group."
        ),
    )
    st.session_state.group_mode = mode

    sample_names = st.session_state.get("sample_names", [])
    st.info(
        f"Found **{len(sample_names)} samples**. "
        "The names below come directly from the loaded table."
    )

    if mode == "Rules by prefix":
        st.markdown("### Rules by prefix")
        st.caption(
            "Generic example: the prefix `A` might correspond to the group `Group 1`. "
            "There are no predefined groups. If there are overlapping prefixes, "
            "the most specific rule is applied first."
        )

        col_sugestao, _ = st.columns([1, 2])
        with col_sugestao:
            if st.button("💡 Suggest based on sample names", width="stretch"):
                sugestao = sugerir_grupos_por_prefixo(sample_names)
                if sugestao:
                    st.session_state.mapping = sugestao
                    st.session_state.pop("mapping_editor", None)
                    st.rerun()
                else:
                    st.info(
                        "No overlapping prefixes found among the samples "
                        "(each name appears unique). Define the groups manually."
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
                    help="Initial part of the sample name."
                ),
                "Group": st.column_config.TextColumn(
                    "Group name", required=False,
                    help="User-defined group name."
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
                    "Detected Group": [preview_group(n) for n in sample_names],
                })
                with st.expander("View how the rules will be applied"):
                    st.dataframe(preview, width="stretch", hide_index=True)
                    missing = int((preview["Detected Group"] == "No group").sum())
                    if missing:
                        st.warning(f"{missing} sample(s) do not correspond to any rule.")

    else:
        st.markdown("### Atribuição das amostras")

        col_sugestao, _ = st.columns([1, 2])
        with col_sugestao:
            if st.button("Suggestion through prefix analysis", width="stretch"):
                prefixos_sugeridos = sugerir_grupos_por_prefixo(sample_names)
                if prefixos_sugeridos:
                    st.session_state.mapping = {
                        nome: grupo
                        for nome in sample_names
                        if (grupo := descobrir_grupo(nome, prefixos_sugeridos)) != "Unknown"
                    }
                    st.session_state.pop("direct_mapping_editor", None)
                    st.rerun()
                else:
                    st.info(
                        "No repeated prefixes found among the samples "
                        "(each name appears to be unique). Assign groups manually."
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
                    "Group", help="Enter the group corresponding to each sample."
                ),
            },
            key="direct_mapping_editor",
        )
        st.caption("You can leave a sample without a group; in that case, it will be marked as `Unknown`.")

    c1, c2 = st.columns(2)
    with c1:
        if st.button("← Return", width="stretch"):
            st.session_state.step = 1
            st.rerun()
    with c2:
        if st.button("Continue →", type="primary", width="stretch"):
            clean = {}
            if mode == "Regras por prefixo":
                for _, row in edited.iterrows():
                    prefix = str(row.get("Prefixo", "")).strip()
                    group = str(row.get("Grupo", "")).strip()
                    if prefix and group and prefix.lower() != "nan" and group.lower() != "nan":
                        clean[prefix] = group
            else:
                for _, row in edited.iterrows():
                    sample = str(row.get("Sample", "")).strip()
                    group = str(row.get("Grupo", "")).strip()
                    if sample and group and group.lower() != "nan":
                        clean[sample] = group

            st.session_state.mapping = clean
            st.session_state.step = 3
            st.rerun()

elif st.session_state.step == 3:
    st.header("3. Review and process")
    st.write("Confirm the data and rules before starting the processing.")

    st.subheader("Files selected")
    for label, key in [
        ("Abundance Table", "abundance"),
        ("Taxonomy Table", "taxonomy"),
        ("FASTA Sequences", "fasta"),
    ]:
        st.write(f"**{label}:** `{st.session_state[key].name}`")

    st.subheader("Configuração dos grupos")
    if st.session_state.mapping:
        st.dataframe(
            pd.DataFrame([
                {"Rule / Sample": k, "Group": v}
                for k, v in st.session_state.mapping.items()
            ]),
            width="stretch",
            hide_index=True,
        )
    else:
        st.warning(
            "No group rules defined. Samples will be processed "
            "and appear as `Unknown` in the metadata."
        )

    st.info(
        "The processing generates the abundance files, TSS abundances, "
        "taxonomy, and metadata."
    )

    c1, c2 = st.columns(2)
    with c1:
        if st.button("← Change Groups", width="stretch"):
            st.session_state.step = 2
            st.rerun()
    with c2:
        if st.button("▶ Execute Processing", type="primary", width="stretch"):
            with st.spinner("Preparing the data…"):
                try:
                    with tempfile.TemporaryDirectory() as tmp:
                        for uploaded in (
                            st.session_state.abundance,
                            st.session_state.taxonomy,
                            st.session_state.fasta,
                        ):
                            uploaded.seek(0)
                        result = process_data(
                            st.session_state.abundance,
                            st.session_state.fasta,
                            st.session_state.taxonomy,
                            st.session_state.mapping,
                            tmp,
                        )
                        result["file_bytes"] = {
                            path.name: Path(path).read_bytes() for path in result["files"]
                        }
                        st.session_state.result = result
                    st.session_state.step = 4
                    st.rerun()
                except Exception as exc:
                    st.error(f"Could not complete the processing:\n\n{exc}")

elif st.session_state.step == 4:
    result = st.session_state.result
    st.header("4. Results")

    if result is None:
        st.warning("There's no result available.")
        if st.button("Back to Start", width="stretch"):
            reset_inputs()
            st.rerun()
        st.stop()

    if result["unknown_samples"]:
        st.warning(
            f"{len(result['unknown_samples'])} sample(s) remained without a group. "
            "They will be identified as `Unknown` in the metadata.txt."
        )
        with st.expander("View samples without a group"):
            st.write(", ".join(result["unknown_samples"]))
    else:
        st.success("✓ All samples were associated with a group.")

    if result["fasta_sem_correspondencia"]:
        with st.expander(
            f"ℹ {len(result['fasta_sem_correspondencia'])} ASV(s) kept the original ID "
            "(No correspondence in FASTA)"
        ):
            st.caption(
                "This is normal if the table already comes with readable IDs instead of sequence hashes. "
                "If you didn't expect this, please confirm that you uploaded the correct FASTA file."
            )
            st.write(", ".join(result["fasta_sem_correspondencia"][:100]))

    if result["fasta_sequencias_duplicadas"]:
        st.warning(
            f"⚠️ {len(result['fasta_sequencias_duplicadas'])} Fasta sequence(s) were"
            "identical to a sequence already processed and were ignored."
        )

    m1, m2, m3 = st.columns(3)
    m1.metric("ASVs", result["asvs"])
    m2.metric("Samples", result["samples"])
    m3.metric("Groups", result["groups"])

    st.subheader("Files generated")
    files_df = pd.DataFrame({
        "File": list(result["file_bytes"].keys()),
        "Description": [
            "Raw abundance table",
            "TSS-normalized table",
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
        file_name="resultados_fase2.zip",
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

    if st.button("↻ Prepare New Data", width="stretch"):
        reset_inputs()
        st.rerun()
