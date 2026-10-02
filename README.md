# MicrobiomeAnalyst Prep


## Graphical interface

This app is a general-purpose data preparation tool for Phase 2. Nothing about your
study is hard-coded: there are no built-in groups, organisms or sample names. Everything
is worked out from the files you upload.

You provide three files:
- an abundance table (`.csv`, `.tsv`, `.xls`, `.xlsx` or `.xlsm`);
- a taxonomy table in any of the same formats;
- a FASTA file.

Once your data is loaded, you can assign samples to groups in one of two ways:
1. **Prefix rules** — a rule applies to every sample whose name starts with the prefix you give it;
2. **Direct assignment** — you pick a group for each sample individually.

Both modes include a **"Suggest from sample names"** button. It proposes groups
automatically, looking only at the sample names in your uploaded file — there are
no predefined groups in the code. The suggestion is always editable before you move on.

You can also choose the normalization method:
- **Relative (TSS)** — the default. Divides each ASV by the sample's total sum.
- **16S gene copy number** — divides each ASV by its most specific taxon's factor,
  using a reference table you upload (e.g. `Taxon` and `CopyNumber` columns).
  Nothing is hard-coded in the app; everything comes from the table you provide.

The app produces four files:
- `otu_table.txt`
- `otu_table_tss.txt` or `otu_table_copynumber.txt` (depending on the chosen method)
- `taxonomy.txt`
- `metadata.txt`

You can download all of them together as a single ZIP.

## Installation

```bash
pip install -r requirements.txt
```

## Running the app

```bash
streamlit run app.py
```

The app will open in your browser.

## Command-line usage

The original command-line workflow still works:

```bash
python scripts/prep.py \
  --input-table data/feature-table.tsv \
  --fasta data/feature.fasta \
  --taxonomy data/taxonomy.tsv \
  --mapping config/mapping.json \
  --output-folder outputs
```

Add `--normalization copy_number --copy-numbers data/copy_numbers.csv` to use
copy-number normalization instead of TSS.

You can keep using `config/mapping.json` to define groups in command-line mode.
The graphical interface doesn't need it.
