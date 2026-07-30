#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: FIN7 
Changes vs FIN6:
1. parse_colabfold_scores() reads ColabFold log.txt and reports
   pLDDT and pTM for the rank_001 (relaxed) model in the run header.
"""

import os, sys, re, subprocess, traceback
import datetime
from pathlib import Path
from collections import Counter

# Konfiguracja Matplotlib bez okna (headless)
os.environ.setdefault("MPLBACKEND", "Agg")

PNG_DPI     = 300
LINE_WRAP   = 60
SEGMENT_LEN = 60
WIDTH_PER_RESIDUE = 0.35

try:
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    mpl.rcParams['savefig.dpi'] = PNG_DPI
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    import logomaker as lm
    from Bio import SeqIO
except ImportError as e:
    print(f"CRITICAL ERROR: Missing Python library: {e}")
    sys.exit(1)

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

LOGO_COLORS = {
    'G': '#00CC00', 'S': '#00CC00', 'T': '#00CC00', 'Y': '#00CC00', 
    'C': '#00CC00', 'Q': '#00CC00', 'N': '#00CC00',
    'K': '#0000CC', 'R': '#0000CC', 'H': '#0000CC',
    'D': '#CC0000', 'E': '#CC0000',
    'A': '#000000', 'V': '#000000', 'L': '#000000', 'I': '#000000', 
    'P': '#000000', 'W': '#000000', 'F': '#000000', 'M': '#000000'
}

# ----------------------- Logger -----------------------

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')
        sys.stdout = self
        sys.stderr = self

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        sys.stdout = self.terminal
        sys.stderr = self.terminal
        self.log.close()

# ----------------------- utils -----------------------

def wrap_iter(s, w):
    for i in range(0, len(s), w):
        yield s[i:i+w]

def wrap_write(fh, s, w):
    for chunk in wrap_iter(s, w):
        fh.write(chunk + "\n")

def parse_range_from_filename(filename):
    m = re.search(r'(\d+)-(\d+)', filename)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)

# ----------------- pipeline steps --------------------

def run_colabfold(input_fasta, output_dir):
    print("=== [1] Running ColabFold ===")
    try:
        subprocess.run(["colabfold_batch", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("ERROR: Command 'colabfold_batch' not found.")
        print("       If you already have results, use the --skip_colabfold flag.")
        raise RuntimeError("ColabFold not found")

    cmd = [
        "colabfold_batch",
        os.path.abspath(input_fasta),
        os.path.abspath(output_dir),
        "--num-recycle", "24",
        "--num-relax", "1"
    ]
    subprocess.run(cmd, check=True)

def run_dssp_in_outdir(output_dir, pdb_basename, dssp_out):
    print("=== [2] Running DSSP via extract_ss_from_pdb.py ===")
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extract_ss_from_pdb.py")
    
    if not os.path.exists(script_path):
        if os.path.exists("extract_ss_from_pdb.py"):
            script_path = "extract_ss_from_pdb.py"
        else:
            raise FileNotFoundError(f"ERROR: 'extract_ss_from_pdb.py' not found.")

    pdb_src = os.path.join(output_dir, pdb_basename)
    if not os.path.exists(pdb_src):
        raise FileNotFoundError(f"PDB not found: {pdb_src}")
    
    cmd = [sys.executable, script_path]
    subprocess.run(cmd, check=True, cwd=output_dir)
    
    produced = os.path.join(output_dir, f"{Path(pdb_basename).stem}.dssp.fasta")
    if not os.path.exists(produced):
        raise FileNotFoundError(f"Expected file not produced: {produced}")
    
    os.replace(produced, dssp_out)
    print(f"[DSSP] -> {dssp_out}")

def get_a3m_file(output_dir):
    candidates = [fn for fn in os.listdir(output_dir) if fn.endswith(".a3m") and "collapsed" not in fn and "filtered" not in fn]
    if not candidates:
        raise FileNotFoundError("No raw .a3m files found in output directory.")
    return os.path.join(output_dir, candidates[0])

def get_query_base_from_a3m(a3m_path):
    return os.path.splitext(os.path.basename(a3m_path))[0]

def filter_a3m_by_coverage(a3m_input, a3m_output, min_cov):
    print(f"=== [FILTER] Filtering sequences shorter than {min_cov*100:.0f}% of query ===")
    
    headers, seqs = [], []
    with open(a3m_input, 'r') as f:
        header = None; buf = []
        for line in f:
            line = line.rstrip("\n")
            if not line: continue
            if line.startswith("#"): continue
            
            if line.startswith(">"):
                if header:
                    seqs.append("".join(buf))
                    buf = []
                header = line
                headers.append(header)
            else:
                buf.append(line)
        if header:
            seqs.append("".join(buf))
            
    if not seqs: raise ValueError("Empty A3M file.")

    query_seq = seqs[0]
    query_len = len(query_seq.replace("-", ""))
    
    kept_count = 0
    total_count = len(seqs)
    
    with open(a3m_output, 'w') as f:
        for h, s in zip(headers, seqs):
            seq_len = len(s.replace("-", ""))
            coverage = seq_len / query_len if query_len > 0 else 0
            
            if h == headers[0] or coverage >= min_cov:
                f.write(h + "\n")
                for chunk in wrap_iter(s, LINE_WRAP):
                    f.write(chunk + "\n")
                kept_count += 1
                
    print(f"[FILTER] Kept {kept_count} of {total_count} sequences (Rejected: {total_count - kept_count})")
    return a3m_output

def collapse_msa(a3m_input, a3m_output):
    print(f"=== [3] Collapsing MSA (Input: {os.path.basename(a3m_input)}) ===")
    headers, seqs = [], []
    with open(a3m_input) as f:
        header=None; buf=[]
        for ln in f:
            ln=ln.rstrip("\n")
            if not ln: continue
            if ln.startswith("#"): continue

            if ln.startswith(">"):
                if header is not None:
                    seqs.append("".join(buf)); buf=[]
                header=ln; headers.append(header)
            else:
                buf.append(ln.replace(" ","").replace("\t",""))
        if header is not None: seqs.append("".join(buf))
    if not seqs: raise ValueError(f"Empty A3M: {a3m_input}")

    def strip_inserts(s): return "".join(ch for ch in s if not ch.islower())
    def sanitize(s):
        out=[]
        for ch in s:
            if ch=='-' or ('A'<=ch<='Z'): out.append(ch)
            else: out.append('-')
        return "".join(out)

    seqs=[sanitize(strip_inserts(s)) for s in seqs]
    query=seqs[0]; L=len(query)

    fixed=[]
    for s in seqs:
        if len(s)<L: s += "-"*(L-len(s))
        elif len(s)>L: s = s[:L]
        fixed.append(s)
    seqs=fixed
    keep=[i for i,a in enumerate(query) if a!='-']
    with open(a3m_output,"w") as f:
        for h,s in zip(headers, seqs):
            collapsed="".join(s[i] for i in keep)
            if h.strip(): f.write(h+"\n")
            wrap_write(f, collapsed, LINE_WRAP)

def convert_to_fasta(a3m_file, fasta_file, query_id=None):
    print("=== [4] Converting A3M to FASTA ===")
    headers, seqs = [], []
    with open(a3m_file) as f:
        header=None; buf=[]
        for ln in f:
            ln=ln.rstrip("\n")
            if not ln: continue
            if ln.startswith("#"): continue

            if ln.startswith(">"):
                if header is not None:
                    seqs.append("".join(buf)); buf=[]
                header=ln; headers.append(header)
            else:
                buf.append(ln)
        if header is not None: seqs.append("".join(buf))
    
    with open(fasta_file,"w") as out:
        for i,(h,s) in enumerate(zip(headers,seqs)):
            hid = query_id if (i==0 and query_id) else (h[1:].split()[0] if h.startswith(">") else f"seq{i}")
            out.write(f">{hid}\n")
            wrap_write(out, s, LINE_WRAP)
    print(f"[fasta] {fasta_file}")

def read_alignment_for_logo(fasta_path):
    seqs=[]
    for rec in SeqIO.parse(str(fasta_path),"fasta"):
        s=str(rec.seq).upper()
        s=re.sub(r"[a-z.]", "", s)
        seqs.append(s)
    if not seqs: raise ValueError("No sequences for logo.")
    Ls={len(s) for s in seqs}
    if len(Ls)!=1: raise ValueError(f"Inconsistent lengths in {fasta_path}: {sorted(Ls)}")
    return seqs

def counts_ignoring_gaps(seqs, alphabet=AA20):
    L=len(seqs[0]); df=pd.DataFrame(0, index=range(L), columns=alphabet, dtype=float)
    for pos in range(L):
        for s in seqs:
            aa=s[pos]
            if aa in alphabet: df.iat[pos, df.columns.get_loc(aa)] += 1.0
    return df

# --------------------- DSSP/SS -----------------------

SS_DESCRIPTIONS = {
    'H': 'α-helix', 'G': '3₁₀-helix', 'I': 'π-helix',
    'E': 'β-strand', 'B': 'β-bridge',
    'T': 'turn', 'S': 'bend', '-': 'coil'
}

def dssp_read_and_fit(ss_fasta, L_expected):
    seqs=[str(r.seq) for r in SeqIO.parse(str(ss_fasta),"fasta")]
    if not seqs: raise ValueError("Empty DSSP FASTA.")
    ss=seqs[0].strip().upper().replace("C","-")
    if len(ss)!=L_expected:
        ss = (ss[:L_expected]) if len(ss)>L_expected else (ss + "-"*(L_expected-len(ss)))
    cnt = Counter(ss)
    print("[DSSP counts]:", " ".join(f"{k}:{cnt[k]}" for k in ['H','G','I','E','B','T','S','-']))
    return ss

def build_palettes(ss_mode="simple", palette="default"):
    if ss_mode == "accurate":
        cats   = ['H','G','I','E','B','T','S','-']
        colors = ['#D95319', '#D95319', '#D95319', '#2878B5', '#8DBDE0', '#9467bd', '#c5b0d5', '#9e9e9e']
        if palette == "reversed":
            colors = ['#C13D0E', '#D95319', '#E37A43', '#0B3B70', '#2878B5', '#9467bd', '#c5b0d5', '#9e9e9e']
        return (cats, colors)
    else:
        cats   = ['Coil','Helix','Sheet']
        if palette == "reversed":
            colors = ['#BFBFBF', '#E37A43', '#4A90D9']
        else:
            colors = ['#BFBFBF', '#D95319', '#2878B5']
        return (cats, colors)

def ss_to_categories(ss_string, mode="simple", palette="default"):
    if mode=="accurate":
        cats, colors = build_palettes("accurate", palette)
        mapping = {c:i for i,c in enumerate(cats)}
        idx = [mapping.get(c, mapping['-']) for c in ss_string]
        legend = [Patch(facecolor=colors[i], edgecolor='none', label=f"{c} {SS_DESCRIPTIONS[c]}") for i,c in enumerate(cats)]
        return np.array(idx, dtype=int), ListedColormap(colors), legend
    else:
        cats, colors = build_palettes("simple", palette)
        def grp(c):
            if c in ('H','G','I'): return 1
            if c in ('E','B'):     return 2
            return 0
        idx = [grp(c) for c in ss_string]
        legend = [
            Patch(facecolor=colors[1], edgecolor='none', label='Helix'),
            Patch(facecolor=colors[2], edgecolor='none', label='Sheet'),
            Patch(facecolor=colors[0], edgecolor='none', label='Coil'),
        ]
        return np.array(idx, dtype=int), ListedColormap(colors), legend

# --------------- rysowanie logo/SS -------------------

def _draw_segment_logo(ax, df_seg, start_idx_abs, tick_step, ylim_bits,
                       start_index_base=1):
    seg_len = len(df_seg)
    logo = lm.Logo(df_seg, ax=ax, 
                   color_scheme=LOGO_COLORS, 
                   font_name='monospace', 
                   font_weight='bold',
                   vpad=0.05,
                   shade_below=0.5, fade_below=0.5)
    
    ax.set_xlim([-0.5, seg_len - 0.5])
    ax.set_ylim(0, ylim_bits)
    ax.set_ylabel("bits", labelpad=5, fontname='monospace', fontsize=10)
    ax.set_xlabel("")
    
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='black', linewidth=1.0)
    
    xt_positions = [i for i in range(seg_len)
                    if tick_step > 0 and ((start_index_base + start_idx_abs + i) % tick_step) == 0]
    xt_labels = [str(start_index_base + start_idx_abs + i) for i in xt_positions]
    ax.set_xticks(xt_positions)
    ax.set_xticklabels(xt_labels, rotation=0, fontname='monospace', fontsize=9)

def _draw_segment_ss(ax_ss, idx_seg, cmap):
    vmax_val = cmap.N - 1
    arr = np.array([idx_seg], dtype=int)
    ax_ss.imshow(arr, aspect='auto', interpolation='none', cmap=cmap,
                 vmin=0, vmax=vmax_val,
                 extent=[-0.5, len(idx_seg)-0.5, 0, 1])
    ax_ss.set_yticks([]); ax_ss.set_xticks([])
    ax_ss.spines[['top','right','left','bottom']].set_visible(False)

def generate_logo_logomaker(aln_fasta, ss_fasta, out_png, title=None,
                            center_last=False, ss_mode="simple", palette="default",
                            start_index_base=1, noise_cutoff=0.0):
    print("=== [5] Generating Logo (Logomaker) ===")
    seqs = read_alignment_for_logo(aln_fasta)
    L = len(seqs[0])
    df_counts = counts_ignoring_gaps(seqs, alphabet=AA20)
    df_info   = lm.transform_matrix(df_counts, from_type='counts', to_type='information')
    
    # === FILTROWANIE SZUMU ===
    if noise_cutoff > 0:
        print(f"[Logo] Filtering noise below {noise_cutoff} bits...")
        df_info[df_info < noise_cutoff] = 0.0
    # =========================
    
    ylim_bits = float(df_info.sum(axis=1).max())
    ylim_bits = max(1.0, round(ylim_bits + 0.1, 2))
    
    ss_string = dssp_read_and_fit(ss_fasta, L_expected=L)
    ss_idx_full, ss_cmap, legend_elems = ss_to_categories(ss_string, mode=ss_mode, palette=palette)

    seg = max(1, SEGMENT_LEN)
    n_rows = (L + seg - 1)//seg
    tick_step = 10 if seg > 60 else 5
    
    fig_w = max(12.0, seg * WIDTH_PER_RESIDUE)
    fig_h = max(6.4, n_rows*3.9 + 2.2)
    
    print(f"[Plotting] Figure Size: {fig_w:.1f} x {fig_h:.1f} inches")
    fig = plt.figure(figsize=(fig_w, fig_h))

    # --- UKŁAD STRONY ---
    if title:
        # Tytuł na samej górze
        fig.suptitle(title, fontsize=24, fontname='monospace', fontweight='bold', y=0.99)

    # Dedykowany pasek legendy: Niżej niż wcześniej (bottom=0.90) i wyższy (height=0.06)
    # [left, bottom, width, height]
    ax_legend_strip = fig.add_axes([0.0, 0.90, 1.0, 0.06])
    ax_legend_strip.axis('off')
    
    # Czcionka powiększona do 18
    leg = ax_legend_strip.legend(handles=legend_elems, ncol=4 if ss_mode=="accurate" else 3,
                                 frameon=False, loc='center', fontsize=18)
    plt.setp(leg.get_texts(), fontname='monospace')

    # Wykresy startują niżej (top=0.85)
    row_logo, row_gap, row_ss = 3.15, 0.5, 0.7
    heights = []
    for _ in range(n_rows):
        heights += [row_logo, row_gap, row_ss, row_gap]

    gs = GridSpec(nrows=len(heights), ncols=1, figure=fig,
                  height_ratios=heights,
                  top=0.85, bottom=0.08, left=0.06, right=0.99, hspace=0.0)

    grid_row = 0
    for r in range(n_rows):
        start = r*seg; end = min((r+1)*seg, L)
        df_seg = df_info.iloc[start:end].reset_index(drop=True)
        ss_idx = ss_idx_full[start:end]; seg_len = end - start
        rem = max(1e-6, seg - seg_len)

        g_logo = gs[grid_row, 0].subgridspec(nrows=1, ncols=2, width_ratios=[seg_len, rem], wspace=0.0)
        ax = fig.add_subplot(g_logo[0, 0])
        
        _draw_segment_logo(ax, df_seg, start_idx_abs=start, tick_step=tick_step,
                           ylim_bits=ylim_bits,
                           start_index_base=start_index_base)
        
        fig.add_subplot(g_logo[0, 1]).axis('off'); grid_row += 2

        g_ss = gs[grid_row, 0].subgridspec(nrows=1, ncols=2, width_ratios=[seg_len, rem], wspace=0.0)
        ax_ss = fig.add_subplot(g_ss[0, 0])
        _draw_segment_ss(ax_ss, ss_idx, ss_cmap)
        fig.add_subplot(g_ss[0, 1]).axis('off'); grid_row += 2

    print(f"Saving to {out_png}...")
    fig.savefig(out_png, dpi=PNG_DPI)
    plt.close(fig)

# -------------------- ColabFold scores ---------------

def parse_colabfold_scores(output_dir):
    """
    Parse ColabFold log.txt and return (pLDDT, pTM) for rank_001 model.
    ColabFold writes a summary line like:
      rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=88.4 pTM=0.82
    Returns (pLDDT_str, pTM_str) or (None, None) if not found.
    """
    log_cf = os.path.join(output_dir, "log.txt")
    if not os.path.exists(log_cf):
        return None, None

    plddt, ptm = None, None
    with open(log_cf, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            # Match summary rank line: "rank_001_... pLDDT=XX.X pTM=0.XX"
            if "rank_001" in line:
                m_plddt = re.search(r'pLDDT=([0-9.]+)', line)
                m_ptm   = re.search(r'pTM=([0-9.]+)', line)
                if m_plddt and m_ptm:
                    plddt = m_plddt.group(1)
                    ptm   = m_ptm.group(1)
                    # Keep the last occurrence (relaxed model line appears after unrelaxed)
    return plddt, ptm


# ------------------------ main -----------------------

def main():
    global SEGMENT_LEN, WIDTH_PER_RESIDUE
    
    if len(sys.argv) < 3:
        print("Usage: python LogoSS_v1.0.py <input_fasta> <output_dir> [options]")
        sys.exit(1)

    input_fasta, output_dir = sys.argv[1], sys.argv[2]
    
    ss_mode, palette = "simple", "default"
    user_start, user_range, start_index_base = None, None, 1
    skip_cf = False
    min_cov = 0.0
    noise_cutoff = 0.05

    extra = sys.argv[3:]
    i = 0
    while i < len(extra):
        a = extra[i]; b = extra[i+1] if (i+1)<len(extra) else None
        if a == "--seg" and b: 
            SEGMENT_LEN = int(b); i += 2
        elif a == "--stretch" and b:
            WIDTH_PER_RESIDUE = float(b); i += 2
        elif a == "--ss" and b in ("simple","accurate"): 
            ss_mode = b; i += 2
        elif a == "--palette" and b in ("default","reversed","cvd"): 
            palette = b; i += 2
        elif a == "--start" and b: 
            user_start = int(b); i += 2
        elif a == "--range" and b and re.match(r"^\d+-\d+$", b):
            user_range = tuple(map(int, b.split("-"))); i += 2
        elif a == "--skip_colabfold":
            skip_cf = True; i += 1
        elif a == "--min_cov" and b:
            min_cov = float(b); i += 2
        elif a == "--noise" and b:
            noise_cutoff = float(b); i += 2
        else: i += 1

    if user_range is not None:
        start_index_base = user_range[0]
    elif user_start is not None:
        start_index_base = user_start
    else:
        fstart, fstop = parse_range_from_filename(os.path.basename(input_fasta))
        if fstart: start_index_base = fstart
    
    os.makedirs(output_dir, exist_ok=True)

    # === [LOGGER START] ===
    log_name = f"{Path(input_fasta).stem}.run.log"
    log_path = os.path.join(output_dir, log_name)
    
    sys_logger = Logger(log_path)
    
    print("="*60)
    print(f"COMMAND: {' '.join(sys.argv)}")
    print(f"CONFIG:  SEG={SEGMENT_LEN}, STRETCH={WIDTH_PER_RESIDUE}, MIN_COV={min_cov}, NOISE={noise_cutoff}")
    print(f"         SS={ss_mode}, PALETTE={palette}, SKIP_CF={skip_cf}")
    print("="*60 + "\n")
    # =======================

    try:
        if not skip_cf:
            run_colabfold(input_fasta, output_dir)
        else:
            print("Skipping ColabFold (--skip_colabfold used)")

        # --- Report ColabFold scores for rank_001 ---
        cf_plddt, cf_ptm = parse_colabfold_scores(output_dir)
        if cf_plddt and cf_ptm:
            print(f"[ColabFold] rank_001 (relaxed)  pLDDT={cf_plddt}  pTM={cf_ptm}")
        else:
            print("[ColabFold] pLDDT/pTM not found in log.txt (scores will be N/A)")
        # -------------------------------------------

        relaxed = sorted([f for f in os.listdir(output_dir)
                          if f.endswith(".pdb") and "_relaxed_rank_001_" in f])
        if not relaxed:
            rank1 = sorted([f for f in os.listdir(output_dir) if f.endswith(".pdb") and "_rank_001_" in f])
            if rank1: best_relaxed = rank1[0]
            else: raise FileNotFoundError("No PDB file found in output directory.")
        else:
            best_relaxed = relaxed[0]

        a3m_input  = get_a3m_file(output_dir)
        query_base = get_query_base_from_a3m(a3m_input)
        
        dssp_out      = os.path.join(output_dir, f"{query_base}.dssp.fasta")
        a3m_collapsed = os.path.join(output_dir, f"{query_base}.collapsed.a3m")
        fasta_output  = os.path.join(output_dir, f"{query_base}.collapsed.fasta")
        logo_png      = os.path.join(output_dir, f"{query_base}.logo.png")

        run_dssp_in_outdir(output_dir, best_relaxed, dssp_out)

        if min_cov > 0:
            a3m_filtered = os.path.join(output_dir, f"{query_base}.filtered.a3m")
            filter_a3m_by_coverage(a3m_input, a3m_filtered, min_cov)
            a3m_input = a3m_filtered

        collapse_msa(a3m_input, a3m_collapsed)
        convert_to_fasta(a3m_collapsed, fasta_output, query_id=query_base)

        logo_title = Path(input_fasta).stem
        
        generate_logo_logomaker(fasta_output, dssp_out, out_png=logo_png,
                                title=logo_title, center_last=False,
                                ss_mode=ss_mode, palette=palette,
                                start_index_base=start_index_base,
                                noise_cutoff=noise_cutoff)

        print("=== SUCCESS ===")
        print(f"Logo generated: {logo_png}")
        print(f"Log saved:      {log_path}")

    except Exception:
        print("\n!!!!!!!!!!!!!! SCRIPT ERROR !!!!!!!!!!!!!!")
        traceback.print_exc()
        sys.exit(1)
    finally:
        sys_logger.close()

if __name__ == "__main__":
    main()
