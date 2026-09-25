import argparse
import json
from pathlib import Path

from prep_core import process_data


def main():
    parser = argparse.ArgumentParser(
        description="Prepara dados do QIIME2 para upload no MicrobiomeAnalyst."
    )
    parser.add_argument("-i", "--tabela-input", required=True,
                        help="Tabela de abundâncias TSV/CSV/Excel.")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Ficheiro FASTA com as sequências ASV.")
    parser.add_argument("-t", "--taxonomia", required=True,
                        help="Ficheiro de taxonomia TSV/CSV/Excel.")
    parser.add_argument("-m", "--mapeamento", required=True,
                        help='JSON com prefixo->grupo. Ex: {"V.CN":"NC","V":"Inoculum"}')
    parser.add_argument("-o", "--pasta-saida", default="output_microbiome",
                        help="Pasta de destino para os ficheiros de output.")
    args = parser.parse_args()

    with open(args.mapeamento, "r", encoding="utf-8") as f:
        dados_grupos = json.load(f)

    result = process_data(
        args.tabela_input,
        args.fasta,
        args.taxonomia,
        dados_grupos,
        Path(args.pasta_saida),
    )

    print("\nFicheiros gerados:")
    for path in result["files"]:
        print(f"   {path}")
    print(f"\nASVs: {result['asvs']}")
    print(f"Amostras: {result['samples']}")
    print(f"Grupos: {result['groups']}")
    if result["unknown_samples"]:
        print(f"AVISO: {len(result['unknown_samples'])} amostra(s) sem grupo: {result['unknown_samples']}")


if __name__ == "__main__":
    main()
