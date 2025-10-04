#!/usr/bin/env python3
"""
signal_extractor_cli_interactive.py

Interactive, user-friendly CLI to extract:
 - signal peptides (N-terminal up to cleavage site)
 - mature proteins (sequence after cleavage)
 - non-signal proteins (full-length)

From a proteome FASTA + SignalP results (prediction_results.txt or output.gff3).

Behavior:
- If tkinter is available and usable, shows a file selection dialog.
- Otherwise, lists found FASTA/GFF/TXT files in current directory for easy selection.
- Allows manual path entry as fallback.
- Streams FASTA; writes outputs incrementally.

Requirements:
    pip install biopython
"""
import os
import glob
import sys
import re
import csv
import argparse

# BioPython
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except Exception as e:
    print("ERROR: Biopython not found. Install with: pip install biopython", file=sys.stderr)
    raise

# Optional tkinter for file dialogs
TK_AVAILABLE = False
try:
    import tkinter as tk
    from tkinter import filedialog, Tk
    TK_AVAILABLE = True
except Exception:
    TK_AVAILABLE = False

# --------------------------
# Utility: file selection UX
# --------------------------
def choose_file_tk(title="Select file", filetypes=(("All files", "*.*"),)):
    """Open a Tk file dialog; return path or None. Safe wrapper."""
    if not TK_AVAILABLE:
        return None
    try:
        root = Tk()
        root.withdraw()
        # try to bring on top (may fail on headless)
        try:
            root.wm_attributes("-topmost", 1)
        except Exception:
            pass
        path = filedialog.askopenfilename(title=title, filetypes=filetypes)
        root.destroy()
        if path:
            return os.path.abspath(path)
        return None
    except Exception:
        return None

def find_files(patterns, search_dir="."):
    """Return sorted list of matching files (non-recursive)."""
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(search_dir, pat)))
    files = sorted(set(files))
    return files

def choose_from_list(files, prompt):
    """Print numbered list and ask user to choose; returns full path or None."""
    if not files:
        return None
    print(prompt)
    for i, f in enumerate(files, start=1):
        print(f"  {i:3d}. {os.path.basename(f)}    ({f})")
    print("  0   Enter a path manually")
    while True:
        choice = input(f"Select file [0-{len(files)}] (press ENTER to retry): ").strip()
        if choice == "":
            continue
        if choice == "0":
            manual = input("Enter full path to file: ").strip()
            if manual and os.path.exists(manual):
                return os.path.abspath(manual)
            print("Path not found. Try again.")
            continue
        try:
            idx = int(choice)
            if 1 <= idx <= len(files):
                return os.path.abspath(files[idx - 1])
        except ValueError:
            pass
        print("Invalid choice. Try again.")

def prompt_for_file(kind="FASTA", extensions=None, prefer_contains=None):
    """
    UX helper: try tkinter -> find files in cwd -> ask manual.
    - kind: label to display
    - extensions: list like ['*.faa','*.fa','*.fasta']
    - prefer_contains: optional substring to prefer in filename (e.g., 'prediction' for SignalP txt)
    """
    print(f"\n--- Select {kind} file ---")
    # try file dialog first
    if TK_AVAILABLE:
        ft = [(f"{kind} files", " ".join(extensions))] if extensions else [("All files", "*.*")]
        path = choose_file_tk(title=f"Select {kind} file", filetypes=ft + [("All files", "*.*")])
        if path:
            print(f"Selected via dialog: {path}")
            return path

    # search cwd
    files = find_files(extensions or ["*.*"], search_dir=".")
    if prefer_contains:
        # move prefered matches to front
        files = sorted(files, key=lambda x: (0 if prefer_contains.lower() in os.path.basename(x).lower() else 1, x))
    if files:
        print(f"Found {len(files)} candidate {kind} files in current directory.")
        sel = choose_from_list(files, f"Choose a {kind} file from the list:")
        if sel:
            return sel

    # manual fallback
    manual = input(f"No {kind} selected. Enter full path manually (or leave blank to cancel): ").strip()
    if manual:
        if os.path.exists(manual):
            return os.path.abspath(manual)
        print("Path not found.")
    return None

# -------------------------
# SignalP parsers (txt/gff)
# -------------------------
def parse_signalp_txt_file(path):
    """Parse SignalP prediction_results.txt -> dict keyed by id_token."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    results = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r'\t+|\s{2,}', line)
        if len(parts) < 2:
            parts = line.split()
            if len(parts) < 2:
                continue
        full_id = parts[0].strip()
        prediction = parts[1].strip()
        id_token = full_id.split()[0]
        cs_left = None
        pr_val = None
        m_cs = re.search(r'CS pos:\s*(\d+)-\d+', line)
        if m_cs:
            cs_left = int(m_cs.group(1))
        m_pr = re.search(r'Pr:\s*([0-9]*\.?[0-9]+)', line)
        if m_pr:
            try:
                pr_val = float(m_pr.group(1))
            except ValueError:
                pr_val = None
        results[id_token] = {
            "full": full_id,
            "prediction": prediction,
            "cs_left": cs_left,
            "pr": pr_val,
            "raw": line,
        }
    return results

def parse_signalp_gff3_file(path):
    """Parse GFF3: return dict seqid -> list of (start,end) coords (1-based inclusive)."""
    mapping = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, feature_type, start, end = cols[0], cols[1], cols[2], cols[3], cols[4]
            if "signal" in feature_type.lower():
                try:
                    s = int(start)
                    e = int(end)
                except Exception:
                    continue
                mapping.setdefault(seqid, []).append((s, e))
    return mapping

# -------------------------
# Extraction (streaming)
# -------------------------
def extract_with_txt(fasta_path, sigp_results, out_prefix, min_pr=0.0, include_pr_in_header=False, summary=False):
    sp_path = f"{out_prefix}_signal_peptides.faa"
    mat_path = f"{out_prefix}_mature_proteins.faa"
    nonsig_path = f"{out_prefix}_nonsignal_proteins.faa"
    summary_path = f"{out_prefix}_summary.csv" if summary else None

    sp_handle = open(sp_path, "w", encoding="utf-8")
    mat_handle = open(mat_path, "w", encoding="utf-8")
    nonsig_handle = open(nonsig_path, "w", encoding="utf-8")
    summary_writer = None
    summary_fh = None
    if summary:
        summary_fh = open(summary_path, "w", newline="", encoding="utf-8")
        summary_writer = csv.writer(summary_fh)
        summary_writer.writerow(["id", "prediction", "cs_left_or_coords", "pr", "original_len", "sp_len", "status"])

    n_sp = n_mat = n_nonsig = 0
    total = 0
    warnings = []

    for rec in SeqIO.parse(fasta_path, "fasta"):
        total += 1
        rec_id = rec.id
        res = sigp_results.get(rec_id)
        # fuzzy match fallback
        if res is None:
            found = None
            for k in sigp_results.keys():
                if k == rec_id or k in rec.description or rec_id in k:
                    found = k
                    break
            if found:
                res = sigp_results[found]
        if res is None:
            SeqIO.write(rec, nonsig_handle, "fasta")
            n_nonsig += 1
            if summary_writer:
                summary_writer.writerow([rec_id, "NO_ENTRY", "", "", len(rec.seq), "", "nonsignal_no_entry"])
            continue
        pred = res.get("prediction", "").upper()
        pr = res.get("pr")
        cs_left = res.get("cs_left")
        if pred == "SP":
            if pr is not None and pr < min_pr:
                SeqIO.write(rec, nonsig_handle, "fasta")
                n_nonsig += 1
                if summary_writer:
                    summary_writer.writerow([rec_id, pred, cs_left or "", pr, len(rec.seq), "", "nonsignal_low_pr"])
                continue
            if cs_left is None:
                warnings.append(f"SP predicted for {rec_id} but no CS pos found; kept as nonsignal.")
                SeqIO.write(rec, nonsig_handle, "fasta")
                n_nonsig += 1
                if summary_writer:
                    summary_writer.writerow([rec_id, pred, "", pr if pr is not None else "", len(rec.seq), "", "nonsignal_no_cs"])
                continue
            if cs_left < 1 or cs_left >= len(rec.seq):
                warnings.append(f"CS pos for {rec_id} out of range (cs_left={cs_left}); kept as nonsignal.")
                SeqIO.write(rec, nonsig_handle, "fasta")
                n_nonsig += 1
                if summary_writer:
                    summary_writer.writerow([rec_id, pred, cs_left, pr if pr is not None else "", len(rec.seq), "", "nonsignal_cs_range"])
                continue
            sp_seq = rec.seq[:cs_left]
            mat_seq = rec.seq[cs_left:]
            desc_sp = f"signal_peptide len={len(sp_seq)}; {res.get('full','')}"
            desc_mat = f"mature_protein len={len(mat_seq)}; {res.get('full','')}"
            if include_pr_in_header and pr is not None:
                desc_sp = f"Pr={pr}; " + desc_sp
                desc_mat = f"Pr={pr}; " + desc_mat
            sp_rec = SeqRecord(Seq(str(sp_seq)), id=rec.id, description=desc_sp)
            mat_rec = SeqRecord(Seq(str(mat_seq)), id=rec.id, description=desc_mat)
            SeqIO.write(sp_rec, sp_handle, "fasta")
            SeqIO.write(mat_rec, mat_handle, "fasta")
            n_sp += 1
            n_mat += 1
            if summary_writer:
                summary_writer.writerow([rec_id, pred, cs_left, pr if pr is not None else "", len(rec.seq), len(sp_seq), "signal"])
        else:
            SeqIO.write(rec, nonsig_handle, "fasta")
            n_nonsig += 1
            if summary_writer:
                summary_writer.writerow([rec_id, pred, "", res.get("pr", ""), len(rec.seq), "", "nonsignal"])

        # small progress print
        if total % 200 == 0:
            print(f"Processed {total} sequences...")

    sp_handle.close()
    mat_handle.close()
    nonsig_handle.close()
    if summary_fh:
        summary_fh.close()

    return {
        "total": total,
        "signal_peptides": n_sp,
        "mature_proteins": n_mat,
        "nonsignal": n_nonsig,
        "warnings": warnings,
        "files": {"sp": sp_path, "mat": mat_path, "nonsig": nonsig_path, "summary": summary_path if summary else None}
    }

def extract_with_gff(fasta_path, gff_map, out_prefix, include_pr_in_header=False, summary=False):
    sp_path = f"{out_prefix}_signal_peptides.faa"
    mat_path = f"{out_prefix}_mature_proteins.faa"
    nonsig_path = f"{out_prefix}_nonsignal_proteins.faa"
    summary_path = f"{out_prefix}_summary.csv" if summary else None

    sp_handle = open(sp_path, "w", encoding="utf-8")
    mat_handle = open(mat_path, "w", encoding="utf-8")
    nonsig_handle = open(nonsig_path, "w", encoding="utf-8")
    summary_writer = None
    summary_fh = None
    if summary:
        summary_fh = open(summary_path, "w", newline="", encoding="utf-8")
        summary_writer = csv.writer(summary_fh)
        summary_writer.writerow(["id", "prediction", "cs_left_or_coords", "pr", "original_len", "sp_len", "status"])

    n_sp = n_mat = n_nonsig = 0
    total = 0
    warnings = []

    for rec in SeqIO.parse(fasta_path, "fasta"):
        total += 1
        seqid = rec.id
        feats = gff_map.get(seqid)
        if feats is None:
            # fuzzy match
            found = None
            for k in gff_map.keys():
                if k == seqid or k in rec.description or seqid in k:
                    found = k
                    break
            if found:
                feats = gff_map[found]
        if not feats:
            SeqIO.write(rec, nonsig_handle, "fasta")
            n_nonsig += 1
            if summary_writer:
                summary_writer.writerow([rec.id, "NO_GFF_FEATURE", "", "", len(rec.seq), "", "nonsignal"])
            continue
        if len(feats) > 1:
            warnings.append(f"{seqid}: multiple signal features in GFF; using first.")
        s, e = feats[0]
        if s < 1 or e > len(rec.seq) or s > e:
            warnings.append(f"{seqid}: GFF coords out of range ({s}-{e}); kept nonsignal.")
            SeqIO.write(rec, nonsig_handle, "fasta")
            n_nonsig += 1
            if summary_writer:
                summary_writer.writerow([rec.id, "GFF_COORD_OOR", f"{s}-{e}", "", len(rec.seq), "", "nonsignal"])
            continue
        sp_seq = rec.seq[s - 1 : e]
        mat_seq = rec.seq[e:]
        desc_sp = f"signal_peptide_gff coords={s}-{e} len={len(sp_seq)}"
        desc_mat = f"mature_protein_gff start={e+1} len={len(mat_seq)}"
        sp_rec = SeqRecord(Seq(str(sp_seq)), id=rec.id, description=desc_sp)
        mat_rec = SeqRecord(Seq(str(mat_seq)), id=rec.id, description=desc_mat)
        SeqIO.write(sp_rec, sp_handle, "fasta")
        SeqIO.write(mat_rec, mat_handle, "fasta")
        n_sp += 1
        n_mat += 1
        if summary_writer:
            summary_writer.writerow([rec.id, "GFF_SIGNAL", f"{s}-{e}", "", len(rec.seq), len(sp_seq), "signal_gff"])
        if total % 200 == 0:
            print(f"Processed {total} sequences...")

    sp_handle.close()
    mat_handle.close()
    nonsig_handle.close()
    if summary_fh:
        summary_fh.close()

    return {
        "total": total,
        "signal_peptides": n_sp,
        "mature_proteins": n_mat,
        "nonsignal": n_nonsig,
        "warnings": warnings,
        "files": {"sp": sp_path, "mat": mat_path, "nonsig": nonsig_path, "summary": summary_path if summary else None}
    }

# -------------------------
# Main interactive flow
# -------------------------
def interactive_flow():
    print("\n=== SignalP Extractor (Interactive) ===\n")
    # select FASTA
    fasta_path = prompt_for_file(kind="proteome FASTA", extensions=["*.faa", "*.fa", "*.fasta"], prefer_contains=None)
    if not fasta_path:
        print("No FASTA selected. Exiting.")
        sys.exit(1)

    # select SignalP result (try txt then gff)
    sigp_path = prompt_for_file(kind="SignalP result (prediction_results.txt or output.gff3)",
                                extensions=["*prediction*.txt", "*.txt", "*.gff3", "*.gff"],
                                prefer_contains="prediction")
    if not sigp_path:
        print("No SignalP result selected. Exiting.")
        sys.exit(1)

    # detect format
    ext = os.path.splitext(sigp_path)[1].lower()
    if ext in (".gff", ".gff3"):
        detected = "gff"
    else:
        # inspect file contents to be safe
        with open(sigp_path, "r", encoding="utf-8", errors="replace") as fh:
            head = "".join([next(fh) for _ in range(30)]) if os.path.getsize(sigp_path) > 0 else ""
        if "CS pos:" in head or head.startswith("# SignalP"):
            detected = "txt"
        else:
            # fallback - if lines with 9 tab cols -> gff
            first_nonblank = ""
            for l in head.splitlines():
                if l.strip():
                    first_nonblank = l
                    break
            if "\t" in first_nonblank and len(first_nonblank.split("\t")) >= 9:
                detected = "gff"
            else:
                detected = "txt"
    print(f"Detected signal file format: {detected.upper()}")

    # options
    while True:
        min_pr_raw = input("Minimum Pr to accept SP (0.0 = accept all) [default 0.0]: ").strip()
        if min_pr_raw == "":
            min_pr = 0.0
            break
        try:
            min_pr = float(min_pr_raw)
            if 0.0 <= min_pr <= 1.0:
                break
        except ValueError:
            pass
        print("Enter a number between 0 and 1 (e.g., 0.9).")

    inc_pr = input("Include Pr value in FASTA headers for SP/mature sequences? (y/N): ").strip().lower() in ("y", "yes")
    out_prefix = input("Output prefix (default 'results'): ").strip()
    if not out_prefix:
        out_prefix = "results"
    summary_opt = input("Write CSV summary file? (Y/n): ").strip().lower()
    summary_flag = False if summary_opt in ("n", "no") else True

    print("\nSummary of choices:")
    print(f"  FASTA : {fasta_path}")
    print(f"  SignalP result : {sigp_path} (format={detected})")
    print(f"  min Pr: {min_pr}")
    print(f"  include Pr in headers: {'Yes' if inc_pr else 'No'}")
    print(f"  out prefix: {out_prefix}")
    print(f"  summary CSV: {'Yes' if summary_flag else 'No'}")

    ok = input("Proceed with extraction? (Y/n): ").strip().lower()
    if ok in ("n", "no"):
        print("Aborted by user.")
        sys.exit(0)

    # parse and extract
    print("\nRunning extraction... this may take a bit for large proteomes.")
    if detected == "txt":
        sigp_results = parse_signalp_txt_file(sigp_path)
        res = extract_with_txt(fasta_path, sigp_results, out_prefix, min_pr=min_pr, include_pr_in_header=inc_pr, summary=summary_flag)
    else:
        gmap = parse_signalp_gff3_file(sigp_path)
        res = extract_with_gff(fasta_path, gmap, out_prefix, include_pr_in_header=inc_pr, summary=summary_flag)

    print("\n--- Extraction completed ---")
    print(f"Total sequences processed: {res['total']}")
    print(f"Signal peptides extracted: {res['signal_peptides']}")
    print(f"Mature proteins generated : {res['mature_proteins']}")
    print(f"Non-signal proteins       : {res['nonsignal']}")
    print("\nFiles written:")
    for k, v in res["files"].items():
        if v:
            print(f"  {k}: {v}")
    if res["warnings"]:
        print("\nWarnings (sample):")
        for w in res["warnings"][:50]:
            print(" - " + w)
    print("\nDone. If you want a GUI instead, run the Streamlit app (streamlit run signal_extractor.py) if you have it.")

# -------------------------
# CLI args support (noninteractive)
# -------------------------
def noninteractive_mode(args):
    # args: parsed Namespace with attributes
    fasta_path = args.fasta
    sigp_path = args.signalp
    if not os.path.exists(fasta_path):
        print(f"ERROR: FASTA not found: {fasta_path}", file=sys.stderr); sys.exit(2)
    if not os.path.exists(sigp_path):
        print(f"ERROR: SignalP file not found: {sigp_path}", file=sys.stderr); sys.exit(2)

    # detect format
    ext = os.path.splitext(sigp_path)[1].lower()
    if ext in (".gff", ".gff3"):
        detected = "gff"
    else:
        with open(sigp_path, "r", encoding="utf-8", errors="replace") as fh:
            head = "".join([next(fh) for _ in range(30)]) if os.path.getsize(sigp_path) > 0 else ""
        if "CS pos:" in head or head.startswith("# SignalP"):
            detected = "txt"
        else:
            first_nonblank = ""
            for l in head.splitlines():
                if l.strip():
                    first_nonblank = l
                    break
            if "\t" in first_nonblank and len(first_nonblank.split("\t")) >= 9:
                detected = "gff"
            else:
                detected = "txt"

    if detected == "txt":
        sigp_results = parse_signalp_txt_file(sigp_path)
        res = extract_with_txt(fasta_path, sigp_results, args.out_prefix, min_pr=args.min_pr, include_pr_in_header=args.include_pr, summary=args.summary)
    else:
        gmap = parse_signalp_gff3_file(sigp_path)
        res = extract_with_gff(fasta_path, gmap, args.out_prefix, include_pr_in_header=args.include_pr, summary=args.summary)

    print("\n--- Extraction summary ---")
    print(f"Total: {res['total']}  SP: {res['signal_peptides']}  NonSP: {res['nonsignal']}")
    for k, v in res["files"].items():
        if v:
            print(f"{k}: {v}")

# -------------------------
# Entrypoint
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="SignalP extractor (interactive-friendly). Run with no args for step-by-step prompts.")
    parser.add_argument("--noninteractive", action="store_true", help="Run in non-interactive (scripted) mode using args below")
    parser.add_argument("--fasta", help="Path to proteome FASTA (required in noninteractive)")
    parser.add_argument("--signalp", help="Path to SignalP result file (txt or gff) (required in noninteractive)")
    parser.add_argument("--min-pr", type=float, default=0.0, help="Minimum Pr to accept SP (0.0 default)")
    parser.add_argument("--include-pr", action="store_true", help="Include Pr in FASTA headers for SP/mature sequences")
    parser.add_argument("--out-prefix", default="results", help="Output prefix (default: results)")
    parser.add_argument("--summary", action="store_true", help="Write CSV summary file")
    args = parser.parse_args()

    if args.noninteractive:
        # require fasta and signalp
        if not args.fasta or not args.signalp:
            parser.error("--noninteractive requires --fasta and --signalp")
        noninteractive_mode(args)
    else:
        interactive_flow()

if __name__ == "__main__":
    main()
