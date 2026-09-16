# LogoSS 0.1.0 (prototype)
Refactored with AI (Claude)

**End-to-end pipeline for annotated sequence logos with secondary structure**

LogoSS takes a single query protein sequence as input and generates a publication-quality sequence logo annotated with secondary structure. It integrates **ColabFold** for homolog retrieval and structure prediction, **DSSP** for secondary-structure annotation, and **Logomaker** for visualization — all in one command.

## What it does

```
Query protein sequence (FASTA)
        ↓
ColabFold (homolog MSA via MMseqs2 + structure prediction)
        ↓
PDB model + a3m alignment + pLDDT/pTM scores
        ↓
DSSP secondary structure          a3m filtering / collapsing
        ↓                                     ↓
        └───────────────┬─────────────────────┘
                        ↓
   Annotated sequence logo (PNG, 300 DPI)
   • per-position information content (Logomaker)
   • α-helix / β-strand / loop track above
```

## Features

- **End-to-end:** single FASTA in → annotated logo PNG out
- **ColabFold integration** — automatic homolog retrieval and structural modeling
- **Skip mode** (`--skip_colabfold`) — re-use existing ColabFold results (in development)
- **MSA filtering by coverage** (`--min_cov`)
- **Two secondary-structure granularities** — simple (H/E/C) or accurate (8 DSSP categories) (in development)
- **Color palette options** — default, reversed, colorblind-friendly (in development)
- **Configurable segmentation** for long proteins (`--seg`, `--stretch`)

## Requirements

- Python 3.9+
- ColabFold (`colabfold_batch` — install via `colabfold-conda` or pip)
- DSSP (`mkdssp` ≥ 4.0)
- Python packages — see `requirements.txt` or `environment.yml`

## Installation

### Option A — conda / mamba / micromamba (recommended)

Installs everything including DSSP in one step:

```bash
git clone https://github.com/<your-username>/LogoSS.git
cd LogoSS
mamba env create -f environment.yml     # or: conda env create -f environment.yml
conda activate logoss
```

ColabFold must be installed separately:
https://github.com/sokrypton/ColabFold#installation

### Option B — pip

```bash
git clone https://github.com/<your-username>/LogoSS.git
cd LogoSS
pip install -r requirements.txt
# Install DSSP and ColabFold separately
```

## Quick usage

```bash
python LogoSS_v1.0.py query.fasta output_dir/
```

That single command will:
1. Run ColabFold on `query.fasta` (writes to `output_dir/`)
2. Run DSSP on the predicted PDB
3. Process the a3m alignment (filter, collapse, convert to FASTA)
4. Generate `logo.png` (300 DPI) with secondary-structure track

### Options

```bash
python LogoSS_v1.0.py query.fasta output_dir/ [options]
```

| Option | Default | Description |
|---|---|---|
| `--skip_colabfold` | off | Re-use existing ColabFold results in `output_dir/` |
| `--seg N` | 60 | Residues per segment in the figure |
| `--stretch X` | 0.35 | Width per residue (inches) |
| `--ss simple\|accurate` | simple | Simple = H/E/C; accurate = 8 DSSP categories (accurate in development) |
| `--palette default\|reversed\|cvd` | default | Color palette (CVD = colorblind-friendly)  (in development) |
| `--min_cov X` | 0.0 | Minimum sequence coverage filter for a3m |
| `--noise_cutoff X` | 0.05 | Letter height cutoff for the logo |

### Example workflow

```bash
# Full pipeline (ColabFold + DSSP + logo)
python LogoSS_v1.0.py PKL_SEAE1.fasta runs/SEAE1/

```

## Citing LogoSS

> ...

LogoSS builds on the following tools — please also cite them:

- **ColabFold** — Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M.
  ColabFold: making protein folding accessible to all. *Nat Methods* 19, 679–682 (2022).
  https://doi.org/10.1038/s41592-022-01488-1
- **AlphaFold2** — Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure
  prediction with AlphaFold. *Nature* 596, 583–589 (2021).
  https://doi.org/10.1038/s41586-021-03819-2
- **MMseqs2** — Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching
  for the analysis of massive data sets. *Nat Biotechnol* 35, 1026–1028 (2017).
  https://doi.org/10.1038/nbt.3988
- **DSSP** — Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition
  of hydrogen-bonded and geometrical features. *Biopolymers* 22, 2577–2637 (1983).
  https://doi.org/10.1002/bip.360221211
- **DSSP 4** — Hekkelman ML, Álvarez Salmoral D, Perrakis A, Joosten RP. DSSP 4: FAIR annotation
  of protein secondary structure. *Protein Sci* 34, e70208 (2025).
  https://doi.org/10.1002/pro.70208
- **Logomaker** — Tareen A, Kinney JB. Logomaker: beautiful sequence logos in Python.
  *Bioinformatics* 36, 2272–2274 (2020). https://doi.org/10.1093/bioinformatics/btz921

## License

LogoSS is released under the **MIT License** — see [LICENSE](LICENSE).

## Author

**Marcin Gradowski**  
Department of Biochemistry and Microbiology  
Warsaw University of Life Sciences — SGGW  
[marcin_gradowski@sggw.edu.pl](mailto:marcin_gradowski@sggw.edu.pl)

## Planned features (v2.0)

The following capabilities are planned for future releases:

- **Custom PDB input** — provide a experimental or pre-computed structure instead of 
  running ColabFold (`--pdb custom.pdb`)
- **Custom secondary structure** — supply a user-defined SS annotation 
  in DSSP-fasta format (`--ss-file custom.dssp.fasta`)
- **Custom alignment** — provide a pre-computed MSA instead of 
  MMseqs2/ColabFold (`--msa custom.fasta` or `.a3m`)
  
Feedback and bug reports are welcome via GitHub Issues.
