#!/usr/bin/env python3
import os, sys, subprocess, tempfile, shutil, glob

CRYST1_LINE = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n".format(
    100.0, 100.0, 100.0, 90.0, 90.0, 90.0
)
# DSSP >=4.0 detects input format from file content: without a HEADER line
# it tries to parse the file as mmCIF and fails with
# "does not seem to be an mmCIF file".
HEADER_LINE = "HEADER    LOGOSS MINIMAL PDB                      01-JAN-00   0000              \n"

def minimal_pdb(src_path, dst_path):
    with open(src_path, "r") as f:
        lines = f.readlines()
    # keep only ATOM/HETATM/TER/END
    atom_block = [ln for ln in lines if ln.startswith(("ATOM", "HETATM", "TER", "END"))]
    if not any(ln.startswith("ATOM") for ln in atom_block):
        raise RuntimeError(f"No ATOM records found in {src_path}")
    with open(dst_path, "w") as g:
        g.write(HEADER_LINE)
        g.write(CRYST1_LINE)
        g.writelines(atom_block)

def run_mkdssp(in_pdb, out_dssp):
    # DSSP >=4.0 dropped the -i/-o short flags in favor of positional
    # input/output arguments, and defaults to mmCIF output unless told
    # otherwise, so --output-format dssp is required for the classic format.
    cmd = ["mkdssp", "--output-format", "dssp", in_pdb, out_dssp]
    subprocess.run(cmd, check=True)

def dssp_to_fasta(dssp_path, fasta_out, header_name):
    mapping = {'H':'H','G':'H','I':'H','E':'E','B':'E','S':'C','T':'C',' ':'C','-':'C'}
    in_table = False
    ss_list = []
    with open(dssp_path, "r") as f:
        for ln in f:
            if not in_table:
                if ln.startswith("  #  RESIDUE"):
                    in_table = True
                continue
            if len(ln) < 17:
                continue
            aa = ln[13]
            ss_char = ln[16]
            if aa in ("!", "*", "X"):
                continue
            ss_list.append(mapping.get(ss_char, 'C'))
    seq = "".join(ss_list)
    wrapped = "\n".join(seq[i:i+70] for i in range(0, len(seq), 70))
    with open(fasta_out, "w") as out:
        out.write(f">{header_name}\n{wrapped}\n")

def process_one(pdb_path):
    base = os.path.basename(pdb_path)
    out_fasta = os.path.splitext(base)[0] + ".dssp.fasta"
    print(f">>> Processing: {base}")
    with tempfile.TemporaryDirectory() as tmp:
        fixed = os.path.join(tmp, "fixed.pdb")
        dssp = os.path.join(tmp, "out.dssp")
        minimal_pdb(pdb_path, fixed)
        run_mkdssp(fixed, dssp)
        dssp_to_fasta(dssp, out_fasta, os.path.splitext(base)[0])
    print(f"[OK] Saved: {out_fasta}")

def main():
    if not shutil.which("mkdssp"):
        print("ERROR: mkdssp not found in PATH", file=sys.stderr)
        sys.exit(2)

    # Use PDB files given as arguments, if any; otherwise take all *.pdb in CWD.
    pdbs = [p for p in sys.argv[1:] if p.endswith(".pdb")]
    if not pdbs:
        pdbs = glob.glob("*.pdb")

    if not pdbs:
        print("No .pdb files found in the working directory.", file=sys.stderr)
        sys.exit(1)

    # Prefer *_relaxed_rank_001_*.pdb files if present; otherwise use all.
    r1 = [p for p in pdbs if "_relaxed_rank_001_" in p]
    targets = r1 if r1 else pdbs

    for p in targets:
        try:
            process_one(p)
        except Exception as e:
            print(f"[ERROR] {p}: {e}", file=sys.stderr)
            # don't abort the whole batch

if __name__ == "__main__":
    main()
