"""Core processing logic for the Bioinformatics Phase 2 application.

This module is deliberately independent from Streamlit so it can be used by
both the graphical interface and the command-line script.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Mapping, Tuple, Union

import pandas as pd

PathLike = Union[str, Path]

# Ordem e prefixos QIIME2/greengenes usados para separar a string de taxonomia.
NIVEIS_TAXONOMICOS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PREFIXOS_TAXONOMICOS = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]


def _open_text(path_or_file):
    """Devolve um stream de texto, tanto para um caminho como para um upload já aberto.

    Não fechamos objetos UploadedFile aqui porque o Streamlit pode precisar
    de os reutilizar dentro da mesma interação.
    """
    if hasattr(path_or_file, "read"):
        return path_or_file, False
    return open(path_or_file, "r", encoding="utf-8"), True


def carregar_mapa_fasta(caminho_fasta) -> Tuple[Dict[str, str], List[str]]:
    """Constrói o mapa md5(sequência) -> ID da ASV a partir de um FASTA.

    Devolve também a lista de IDs cuja sequência colide com a de um ID
    anterior (mesma sequência, cabeçalho diferente), para o utilizador poder
    ser avisado em vez de a colisão passar despercebida.
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
        raise ValueError("O ficheiro FASTA não contém sequências válidas.")
    return mapa, colisoes


def limpar_taxonomia(texto_tax):
    """Separa uma string tipo 'k__Bacteria;p__Firmicutes;...' pelos seus níveis."""
    resultado = {nivel: "" for nivel in NIVEIS_TAXONOMICOS}
    if pd.isna(texto_tax):
        return resultado
    for parte in str(texto_tax).split(";"):
        parte = parte.strip()
        for prefixo, nivel in zip(PREFIXOS_TAXONOMICOS, NIVEIS_TAXONOMICOS):
            if parte.startswith(prefixo):
                resultado[nivel] = parte[len(prefixo):].strip()
    return resultado


def formatar_taxonomia(df_tax: pd.DataFrame) -> pd.DataFrame:
    """Devolve uma tabela Kingdom..Species a partir de dois formatos possíveis:

    - uma única coluna com strings 'k__x;p__y;...' (export QIIME2/greengenes); ou
    - colunas já separadas com os nomes dos níveis (ex.: export do MicrobiomeAnalyst).

    Se nenhuma coluna corresponder a um nível conhecido, assume-se o primeiro
    formato e usa-se a primeira coluna da tabela.
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
    """Determina o grupo pelo prefixo mais longo que corresponde (ex.: V.CN antes de V)."""
    for prefixo in sorted(dicionario_grupos.keys(), key=len, reverse=True):
        if nome_amostra.startswith(prefixo):
            return dicionario_grupos[prefixo]
    return "Unknown"


def _detetar_separador(primeira_linha: str) -> str:
    """Tenta perceber se a linha usa vírgulas ou tabs, contando ocorrências."""
    if primeira_linha.count(",") > primeira_linha.count("\t"):
        return ","
    return "\t"


def _ler_primeiras_linhas(file_or_path, n: int = 2) -> List[str]:
    """Lê as primeiras `n` linhas de um ficheiro ou upload, sem consumir o stream."""
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
    """Lê CSV/TSV/XLS/XLSX, detetando o separador pelo conteúdo em vez de assumir pela extensão.

    Também reconhece o cabeçalho de comentário que o QIIME2 coloca antes da
    linha de colunas nas tabelas de abundância exportadas.
    """
    nome = getattr(file_or_path, "name", str(file_or_path)).lower()
    sufixo = Path(nome).suffix.lower()

    if sufixo in {".xls", ".xlsx", ".xlsm"}:
        return pd.read_excel(file_or_path, index_col=0)

    linhas = _ler_primeiras_linhas(file_or_path)
    primeira_linha = linhas[0] if linhas else ""

    separador = "," if sufixo == ".csv" else _detetar_separador(primeira_linha)

    # Ficheiros exportados do QIIME2 têm uma linha "# Constructed from biom
    # file" antes do cabeçalho real. Só ignoramos essa linha se não parecer
    # já ser o próprio cabeçalho (que também começa por '#').
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
    """Valida os requisitos estruturais mínimos antes de processar os dados."""
    errors = []
    if df_otu.empty:
        errors.append("A tabela de abundâncias está vazia.")
    if df_tax.empty:
        errors.append("A tabela de taxonomia está vazia.")
    if df_otu.index.duplicated().any():
        errors.append("A tabela de abundâncias contém IDs de ASV duplicados.")
    if df_tax.index.duplicated().any():
        errors.append("A tabela de taxonomia contém IDs de ASV duplicados.")
    if df_otu.shape[1] == 0:
        errors.append("Não foram encontradas amostras na tabela de abundâncias.")

    try:
        numeric = df_otu.apply(pd.to_numeric)
        if numeric.isna().any().any():
            errors.append("A tabela de abundâncias contém valores que não são numéricos.")
    except Exception:
        errors.append("Não foi possível interpretar a tabela de abundâncias como valores numéricos.")

    if not fasta_map:
        errors.append("Não foram encontradas sequências no FASTA.")
    return errors


def process_data(
    tabela_input,
    fasta,
    taxonomia,
    mapeamento: Mapping[str, str],
    pasta_saida: PathLike,
) -> dict:
    """Corre o pipeline completo da Fase 2 e devolve um resumo do resultado."""
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

    # Quantas ASVs da tabela não tinham correspondência no FASTA? Se o índice
    # já vier com os IDs em vez de hashes de sequência, é normal que seja
    # tudo — mas vale a pena reportar para o utilizador perceber o que aconteceu.
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
        raise ValueError(f"Existem amostras sem contagens (soma igual a zero): {bad}")

    otu_tss = df_otu.div(col_sums, axis=1)
    desvio_max = float((otu_tss.sum(axis=0) - 1.0).abs().max())
    if desvio_max >= 1e-9:
        raise ValueError(f"Validação TSS falhou (desvio máximo: {desvio_max:.2e}).")
    otu_tss.index.name = "#NAME"
    tss_path = output_dir / "otu_table_tss.txt"
    otu_tss.to_csv(tss_path, sep="\t", float_format="%.15g")

    if df_tax.shape[1] == 0:
        raise ValueError("A tabela de taxonomia não contém nenhuma coluna de taxonomia.")
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
