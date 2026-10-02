import argparse
import json
from pathlib import Path

from prep_core import process_data


def main():
    parser = argparse.ArgumentParser(
        description="Prepares QIIME2 data for upload to MicrobiomeAnalyst."
    )
    parser.add_argument("-i", "--input-table", required=True,
                        help="Abundance table, TSV/CSV/Excel.")
    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA file with the ASV sequences.")
    parser.add_argument("-t", "--taxonomy", required=True,
                        help="Taxonomy file, TSV/CSV/Excel.")
    parser.add_argument("-m", "--mapping", required=True,
                        help='JSON with prefix->group. Ex: {"A.CN":"Group1","A":"Group2"}')
    parser.add_argument("-o", "--output-folder", default="output_microbiome",
                        help="Destination folder for the output files.")
    parser.add_argument("-n", "--normalization", choices=["tss", "copy_number"], default="tss",
                        help="Normalization method (default: tss).")
    parser.add_argument("-c", "--copy-numbers",
                        help="Copy number reference table (required if --normalization copy_number).")
    args = parser.parse_args()

    with open(args.mapping, "r", encoding="utf-8") as f:
        group_mapping = json.load(f)

    result = process_data(
        args.input_table,
        args.fasta,
        args.taxonomy,
        group_mapping,
        Path(args.output_folder),
        normalization_method=args.normalization,
        copy_number_table=args.copy_numbers,
    )

    print("\nGenerated files:")
    for path in result["files"]:
        print(f"   {path}")
    print(f"\nASVs: {result['asvs']}")
    print(f"Samples: {result['samples']}")
    print(f"Groups: {result['groups']}")
    if result["unknown_samples"]:
        print(f"WARNING: {len(result['unknown_samples'])} sample(s) without a group: {result['unknown_samples']}")


if __name__ == "__main__":
    main()
