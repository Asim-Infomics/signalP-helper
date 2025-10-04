**SignalP-Helper (CLI + Streamlit)**
SignalP-Helper is a lightweight, user-friendly tool designed for biologists to extract signal peptides, mature proteins, and non-signal proteins from large proteome FASTA files using SignalP 6.x prediction results.

Instead of manually cross-referencing SignalP output files with your FASTA sequences (a time-consuming and error-prone task), this tool automates the process in a simple workflow.

**Input**

Proteome FASTA file (.faa, .fa, .fasta)

SignalP prediction results (prediction_results.txt or output.gff3)

**Output**

*_signal_peptides.faa — N-terminal signal peptide sequences (up to cleavage site)

*_mature_proteins.faa — sequences after cleavage site

*_nonsignal_proteins.faa — proteins without predicted signals

*_summary.csv (optional) — tabular summary of predictions, cleavage sites, and probabilities

**Features**

Interactive CLI: step-by-step file selection, no need for complex commands.

Beginner-friendly: designed for biologists unfamiliar with directories or scripting.

Flexible parsing: supports both prediction_results.txt and output.gff3 formats.

Efficient: streams FASTA input without loading everything into memory.

Biologist-friendly outputs: FASTA files separated into signal, mature, and nonsignal proteins with optional CSV summary.

Optional GUI: a Streamlit version for researchers preferring a graphical interface.

**Installation**
Create a Python environment (recommended):

python -m venv venv
source venv/bin/activate    # Linux/macOS
# venv\Scripts\activate     # Windows
pip install -r requirements.txt


requirements.txt contains:

biopython
streamlit   # optional, for GUI
pandas      # optional, for summary CSV
matplotlib  # optional, for visualization


If you only need the CLI, biopython is sufficient.

**Usage**

Interactive mode (recommended for most users):
python signal_extractor_cli_interactive.py


This will guide you step by step to select your FASTA and SignalP result files.

Non-interactive mode (for advanced users):
python signal_extractor_cli_interactive.py --noninteractive \
  --fasta proteome.faa \
  --signalp prediction_results.txt \
  --min-pr 0.9 \
  --out-prefix myresults \
  --summary

Example outputs

If you set --out-prefix myresults:

myresults_signal_peptides.faa

myresults_mature_proteins.faa

myresults_nonsignal_proteins.faa

myresults_summary.csv (if --summary is enabled)

Streamlit GUI (optional)

If you prefer a graphical interface:

streamlit run signal_extractor.py


This allows file uploads, inspection of warnings, and direct downloads.

Notes & Troubleshooting

Cleavage site parsing: SignalP reports cleavage sites like CS pos: 22-23. The tool interprets the left position (22) as the signal peptide length.

FASTA header matching: Matching is done by the record ID (first whitespace-separated token). Fuzzy matching is applied if IDs differ slightly.

Uncertain predictions: If a signal peptide is predicted but has no valid cleavage site, the sequence is placed into nonsignal_proteins.faa and logged as a warning.

**Why this tool?**
Biologists often spend hours manually matching SignalP predictions with raw sequences. This project was created during MPhil research on Puccinia striiformis to save time and provide a reliable workflow for researchers handling large proteomes.
