# design-filter

Filter computational protein designs by hotspot–binder atomic contacts.

Given a set of PDB files (e.g. from RFdiffusion + ProteinMPNN + AlphaFold2),
this tool identifies designs where specified **hotspot residues** on the target
chain make close contacts with atoms on the binder chain(s). Designs that pass
the distance cutoff are copied to an output directory, and a CSV summary of all
contacts is generated.

## Quick start

```bash
pip install -r requirements.txt

# 1. Put your PDB files in a directory (default: designs/)
mkdir designs
cp /path/to/your/*.pdb designs/

# 2. Edit config.txt to match your target
#    chain    – target chain ID
#    residues – hotspot residue numbers (comma-separated)
#    cutoff   – distance threshold in angstroms

# 3. Run
python filter.py
```

Passing designs are copied to `designs_filtered/` and a contact summary is
written to `summary.csv`.

## Try with example data

```bash
python filter.py -i examples
```

## Config file format

```
# Target chain ID
chain=A

# Hotspot residue numbers (comma-separated)
residues=86,90,94

# Distance cutoff in angstroms
cutoff=4.5
```

## Command-line options

```
python filter.py [OPTIONS]

  -c, --config      Path to config file          (default: config.txt)
  -i, --input-dir   Directory of input PDB files  (default: designs/)
  -o, --output-dir  Directory for passing designs  (default: designs_filtered/)
  -s, --summary     Output CSV path               (default: summary.csv)
```

## Output

`summary.csv` contains one row per passing design:

| filename | contact_count | contact_1 | distance_1_A | contact_2 | distance_2_A | ... |
|----------|---------------|-----------|--------------|-----------|--------------|-----|
| design_1.pdb | 3 | A.86.CA - B.12.CB | 3.72 | A.86.N - B.12.O | 4.01 | ... |

## Requirements

- Python 3.8+
- BioPython

## License

MIT
