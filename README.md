# LogoSS

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
        ↓                                    ↓
        └───────────────┬─────────────────────┘
                        ↓
   Annotated sequence logo (PNG, 300 DPI)
   • per-position information content (Logomaker)
   • α-helix / β-strand / loop track above
   • pLDDT / pTM reported in figure header
```

## Features

- **End-to-end:** single FASTA in → annotated logo PNG out
- **ColabFold integration** — automatic homolog retrieval and structural modeling
- **Skip mode** (`--skip_colabfold`) — re-use existing ColabFold results
- **MSA filtering by coverage** (`--min_cov`)
- **Two secondary-structure granularities** — simple (H/E/C) or accurate (8 DSSP categories)
- **Color palette options** — default, reversed, colorblind-friendly
- **Configurable segmentation** for long proteins (`--seg`, `--stretch`)
- Reports **pLDDT and pTM** of the representative model in the figure title

## Requirements

- Python 3.9+
- ColabFold (`colabfold_batch` — install via `colabfold-conda` or pip)
- DSSP (`mkdssp` ≥ 3.0)
- Python packages — see `requirements.txt` or `environment.yml`

## Installation

### Option A — conda / mamba (recommended)

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
| `--ss simple\|accurate` | simple | Simple = H/E/C; accurate = 8 DSSP categories |
| `--palette default\|reversed\|cvd` | default | Color palette (CVD = colorblind-friendly) |
| `--min_cov X` | 0.0 | Minimum sequence coverage filter for a3m |
| `--noise_cutoff X` | 0.05 | Letter height cutoff for the logo |

### Example workflow

```bash
# Full pipeline (ColabFold + DSSP + logo)
python LogoSS_v1.0.py PKL_SEAE1.fasta runs/SEAE1/

# Re-render figure with different segmentation (no re-prediction)
python LogoSS_v1.0.py PKL_SEAE1.fasta runs/SEAE1/ \
    --skip_colabfold --seg 80 --palette cvd
```

## Citing LogoSS

If you use LogoSS, please cite:

> Gradowski M. (2026). *Systematic discovery of protein kinase-like domains reveals diverse evolutionary strategies in the human oral microbiome.* [Journal, DOI]

A `CITATION.cff` is included for direct GitHub citation export.

## License

LogoSS is released under the **MIT License** — see [LICENSE](LICENSE).

## Author

**Marcin Gradowski**  
Department of Biochemistry and Microbiology  
Warsaw University of Life Sciences — SGGW  
[marcin_gradowski@sggw.edu.pl](mailto:marcin_gradowski@sggw.edu.pl)

## Status

LogoSS is at **version 1.0**. Feedback and bug reports are welcome via GitHub Issues.
