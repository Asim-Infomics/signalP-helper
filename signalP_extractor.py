# signal_extractor.py
"""
Streamlit app to extract signal peptides and mature proteins from
a FASTA proteome using SignalP output (prediction_results.txt or output.gff3).

Usage:
    pip install streamlit biopython
    streamlit run signal_extractor.py
"""

import re
from io import BytesIO, StringIO
from typing import Dict, List, Optional, Tuple

import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_signalp_txt(text: str) -> Dict[str, Dict]:
    """
    Parse SignalP prediction_results.txt style output.
    Returns a dict keyed by protein ID_token (first token), value with:
      { 'full_id_line': str,
        'prediction': 'SP' or 'OTHER' or ...,
        'cs_left': int or None,  # cleavage occurs between cs_left and cs_left+1
        'pr': float or None,
        'raw': original line
      }
    """
    results = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # split by tabs OR run of 2+ spaces
        parts = re.split(r'\t+|\s{2,}', line)
        if len(parts) < 2:
            continue
        full_id = parts[0].strip()
        prediction = parts[1].strip()
        # compute id_token as the first whitespace-separated token (typical FASTA id)
        id_token = full_id.split()[0]
        cs_left = None
        pr_val = None
        # find CS pos pattern like: "CS pos: 22-23" or "CS pos: 22-23. Pr: 0.9811"
        m_cs = re.search(r'CS pos:\s*(\d+)-\d+', line)
        if m_cs:
            cs_left = int(m_cs.group(1))
        m_pr = re.search(r'Pr:\s*([0-9]*\.?[0-9]+)', line)
        if m_pr:
            try:
                pr_val = float(m_pr.group(1))
            except:
                pr_val = None
        results[id_token] = {
            'full_id_line': full_id,
            'prediction': prediction,
            'cs_left': cs_left,
            'pr': pr_val,
            'raw': line
        }
    return results


def parse_signalp_gff3(text: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse GFF3 and return mapping of seqid -> list of (start, end) for signal peptide features.
    Looks for feature type containing 'signal' (case-insensitive).
    GFF coords are 1-based inclusive.
    """
    mapping = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        seqid, source, feature_type, start, end, score, strand, phase, attributes = cols[:9]
        if 'signal' in feature_type.lower():
            try:
                s = int(start)
                e = int(end)
            except:
                continue
            mapping.setdefault(seqid, []).append((s, e))
    return mapping


def extract_sequences(
    records: List[SeqRecord],
    sigp_results: Dict[str, Dict],
    min_pr: float = 0.0
) -> Tuple[List[SeqRecord], List[SeqRecord], List[SeqRecord], List[str]]:
    """
    Using parsed txt results, extract:
      - signal peptide records (id and seq)
      - mature protein records
      - nonsignal (full-length) protein records
    Returns lists and a list of warnings (missing ids etc).
    cleavage: uses cs_left as number of residues in signal peptide (so SP = seq[:cs_left])
    """
    sp_records = []
    mature_records = []
    nonsignal_records = []
    warnings = []

    # build dict of fasta records by id_token (first token)
    fasta_by_id = {rec.id: rec for rec in records}

    for rec in records:
        rec_id = rec.id
        res = sigp_results.get(rec_id)
        if not res:
            # try fuzzy match: find any result key that is substring of rec.description or vice versa
            matches = [k for k, v in sigp_results.items() if k in rec.description or k == rec.id]
            if matches:
                res = sigp_results[matches[0]]
        if not res:
            # No prediction entry -> treat as nonsignal but warn
            nonsignal_records.append(rec)
            warnings.append(f"No SignalP entry found for {rec_id}; kept as nonsignal.")
            continue

        pred = res.get('prediction', '').upper()
        if pred == 'SP':
            pr = res.get('pr')
            if pr is not None and pr < min_pr:
                # treat as nonsignal due to low confidence
                nonsignal_records.append(rec)
                continue
            cs_left = res.get('cs_left')
            if cs_left is None:
                warnings.append(f"SP predicted for {rec_id} but no CS pos found; kept as nonsignal.")
                nonsignal_records.append(rec)
                continue
            if cs_left < 1 or cs_left >= len(rec.seq):
                warnings.append(f"CS pos for {rec_id} seems out of range (cs_left={cs_left}). Kept as nonsignal.")
                nonsignal_records.append(rec)
                continue
            # signal peptide: residues 0..cs_left-1 ; mature: cs_left..end
            sp_seq = rec.seq[:cs_left]
            mat_seq = rec.seq[cs_left:]
            sp_rec = SeqRecord(sp_seq, id=rec.id, description=f"signal_peptide len={len(sp_seq)}; {res.get('full_id_line','')}")
            mat_rec = SeqRecord(mat_seq, id=rec.id, description=f"mature_protein len={len(mat_seq)}; {res.get('full_id_line','')}")
            sp_records.append(sp_rec)
            mature_records.append(mat_rec)
        else:
            nonsignal_records.append(rec)

    return sp_records, mature_records, nonsignal_records, warnings


def extract_from_gff3(
    records: List[SeqRecord],
    gff_map: Dict[str, List[Tuple[int, int]]]
) -> Tuple[List[SeqRecord], List[SeqRecord], List[SeqRecord], List[str]]:
    """
    Use GFF3 mapping seqid->[(start,end)] to extract SPs.
    If there are multiple SP features per seq, we use the first one (warn).
    """
    sp_records = []
    mature_records = []
    nonsignal_records = []
    warnings = []

    fasta_dict = {rec.id: rec for rec in records}
    for rec in records:
        seqid = rec.id
        feats = gff_map.get(seqid, [])
        if not feats:
            nonsignal_records.append(rec)
            continue
        if len(feats) > 1:
            warnings.append(f"{seqid}: multiple signal features in GFF3; using the first one.")
        s, e = feats[0]
        # GFF is 1-based inclusive; signal peptide segment = seq[s-1:e]
        if s < 1 or e > len(rec.seq) or s > e:
            warnings.append(f"{seqid}: GFF coordinates out of range ({s}-{e}); kept as nonsignal.")
            nonsignal_records.append(rec)
            continue
        sp_seq = rec.seq[s - 1: e]
        # assume cleavage after end, so mature starts at e+1 (index e)
        mature_seq = rec.seq[e:]
        sp_rec = SeqRecord(sp_seq, id=rec.id, description=f"signal_peptide_gff len={len(sp_seq)}; coords={s}-{e}")
        mat_rec = SeqRecord(mature_seq, id=rec.id, description=f"mature_protein_gff len={len(mature_seq)}; coords={e+1}-{len(rec.seq)}")
        sp_records.append(sp_rec)
        mature_records.append(mat_rec)
    return sp_records, mature_records, nonsignal_records, warnings


def fasta_to_bytes(records: List[SeqRecord]) -> bytes:
    bio = StringIO()
    SeqIO.write(records, bio, "fasta")
    return bio.getvalue().encode()


# --- Streamlit UI ---
st.set_page_config(page_title="SignalP extractor", layout="wide")
st.title("SignalP → Signal peptide & mature protein extractor")
st.markdown(
    """
Upload a FASTA proteome (e.g. `.faa`) and a SignalP **prediction_results.txt** file or a **GFF3** output.
This tool will extract:
- signal peptides (N-terminal up to cleavage site)
- mature proteins (sequence after cleavage)
- non-signal proteins (full-length sequences for proteins not predicted as SP)

**Notes**
- For SignalP `CS pos: X-Y` the app uses **X** as the signal peptide length (cleavage between X and X+1).
- If you have `Pr:` values, you can filter by minimum confidence.
"""
)

col1, col2 = st.columns(2)

with col1:
    fasta_file = st.file_uploader("Upload FASTA file (proteome) (.faa/.fa/.fasta)", type=["faa", "fa", "fasta"], help="Example: .faa")
    sigp_file = st.file_uploader("Upload SignalP result (prediction_results.txt or output.gff3)", type=["txt", "gff3", "gff"], help="Prediction results from SignalP")

with col2:
    st.write("Options")
    min_pr = st.number_input("Min Pr (confidence) to accept SP (0.0 = accept all) — if present in file", min_value=0.0, max_value=1.0, value=0.0, step=0.01)
    show_warnings = st.checkbox("Show parsing warnings", value=True)
    parse_method = st.radio("Interpret signal file as:", ("Auto-detect", "SignalP TXT (prediction_results.txt)", "GFF3"))

process = st.button("Extract sequences")

if process:
    if not fasta_file:
        st.error("Please upload a FASTA file first.")
    elif not sigp_file:
        st.error("Please upload a SignalP result file.")
    else:
        # read fasta
        try:
            fasta_text = fasta_file.read().decode()
        except Exception:
            # if binary we still try .read and use BytesIO
            fasta_file.seek(0)
            fasta_text = fasta_file.read().decode('latin1')
        # parse fasta into SeqRecord list
        try:
            records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
            if not records:
                st.error("No sequences found in the FASTA file. Check format.")
                st.stop()
        except Exception as e:
            st.error(f"Failed to parse FASTA: {e}")
            st.stop()

        # read signalp file
        try:
            sigp_text = sigp_file.read().decode()
        except Exception:
            sigp_file.seek(0)
            sigp_text = sigp_file.read().decode('latin1')

        # determine parsing approach
        use_gff = False
        if parse_method == "GFF3":
            use_gff = True
        elif parse_method == "SignalP TXT (prediction_results.txt)":
            use_gff = False
        else:
            # Auto-detect
            first_nonblank = ""
            for s in sigp_text.splitlines():
                s = s.strip()
                if s:
                    first_nonblank = s
                    break
            if first_nonblank.startswith("# SignalP") or "CS pos:" in sigp_text or re.search(r'\t', first_nonblank):
                use_gff = False
            else:
                # simple heuristic: GFF3 lines are tab-delimited with 9 columns and not starting with '#'
                cols = first_nonblank.split("\t")
                if len(cols) >= 9 and cols[2].strip():
                    use_gff = True
                else:
                    # default to txt
                    use_gff = False

        st.info(f"Parsing as {'GFF3' if use_gff else 'SignalP TXT'} results")

        if use_gff:
            gmap = parse_signalp_gff3(sigp_text)
            sp_recs, mat_recs, nonsig_recs, warnings = extract_from_gff3(records, gmap)
        else:
            sigp_results = parse_signalp_txt(sigp_text)
            sp_recs, mat_recs, nonsig_recs, warnings = extract_sequences(records, sigp_results, min_pr=min_pr)

        # results summary
        st.subheader("Results summary")
        st.write(f"Total sequences in FASTA: **{len(records)}**")
        st.write(f"Signal peptides extracted: **{len(sp_recs)}**")
        st.write(f"Mature proteins (from SP): **{len(mat_recs)}**")
        st.write(f"Non-signal proteins (full-length): **{len(nonsig_recs)}**")

        if show_warnings and warnings:
            st.warning("Warnings (parsing & edge cases):")
            st.write("\n".join(warnings[:50]))
            if len(warnings) > 50:
                st.write(f"... and {len(warnings)-50} more")

        # show sample sequences
        def show_sample(title, recs, n=5):
            st.write(f"**{title}** (showing up to {n}):")
            for r in recs[:n]:
                seq_preview = str(r.seq)[:200]
                st.code(f">{r.id} {r.description}\n{seq_preview}{'...' if len(r.seq)>200 else ''}")
            if not recs:
                st.write("_none_")

        st.markdown("---")
        sample_cols = st.columns(3)
        with sample_cols[0]:
            show_sample("Signal peptides", sp_recs)
        with sample_cols[1]:
            show_sample("Mature proteins", mat_recs)
        with sample_cols[2]:
            show_sample("Non-signal proteins", nonsig_recs)

        # provide downloads
        st.markdown("---")
        st.subheader("Download FASTA files")
        if sp_recs:
            sp_bytes = fasta_to_bytes(sp_recs)
            st.download_button("Download signal_peptides.faa", sp_bytes, file_name="signal_peptides.faa", mime="text/plain")
        else:
            st.info("No signal peptides to download.")

        if mat_recs:
            mat_bytes = fasta_to_bytes(mat_recs)
            st.download_button("Download mature_proteins.faa", mat_bytes, file_name="mature_proteins.faa", mime="text/plain")
        else:
            st.info("No mature proteins to download.")

        if nonsig_recs:
            non_bytes = fasta_to_bytes(nonsig_recs)
            st.download_button("Download nonsignal_proteins.faa", non_bytes, file_name="nonsignal_proteins.faa", mime="text/plain")
        else:
            st.info("No nonsignal proteins to download.")

        st.success("Done — files ready for download.")
