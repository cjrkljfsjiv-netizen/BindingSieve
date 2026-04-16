#!/usr/bin/env python3
"""
Filter protein designs by hotspot-binder atomic interactions.

Scans PDB files in an input directory, identifies designs where specified
hotspot residues on a target chain make close contacts (within a distance
cutoff) with atoms on other (binder) chains, copies passing designs to an
output directory, and writes a CSV summary of all contacts found.

Usage:
    python filter.py                        # uses default config.txt
    python filter.py -c my_config.txt       # custom config file
    python filter.py -i designs -o passed   # override directories
"""

import argparse
import csv
import os
import shutil
import sys

from Bio.PDB import PDBParser


# ── Config ───────────────────────────────────────────────────────────────

def load_config(path):
    """Read a key=value config file and return a dict."""
    cfg = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                key, value = line.split("=", 1)
                cfg[key.strip()] = value.strip()
    return cfg


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", default="config.txt",
                    help="Path to config file (default: config.txt)")
    ap.add_argument("-i", "--input-dir", default=None,
                    help="Directory of input PDB files (default: designs/)")
    ap.add_argument("-o", "--output-dir", default=None,
                    help="Directory for passing designs (default: designs_filtered/)")
    ap.add_argument("-s", "--summary", default="summary.csv",
                    help="Output CSV path (default: summary.csv)")
    return ap.parse_args()


# ── Analysis ─────────────────────────────────────────────────────────────

def find_contacts(structure, target_chain_id, hotspot_resnums, cutoff):
    """
    Return a list of (contact_label, distance) tuples for atoms on the
    target chain's hotspot residues that are within *cutoff* angstroms of
    any atom on another chain.
    """
    model = structure[0]

    try:
        target_chain = model[target_chain_id]
    except KeyError:
        return []

    # Collect hotspot atoms
    hotspot_atoms = []
    for res in target_chain:
        if res.id[1] in hotspot_resnums:
            hotspot_atoms.extend((res, atom) for atom in res)
    if not hotspot_atoms:
        return []

    # Collect binder atoms (all chains except target)
    binder_atoms = []
    for chain in model:
        if chain.id == target_chain_id:
            continue
        for res in chain:
            binder_atoms.extend((chain, res, atom) for atom in res)

    # Pairwise distance check
    contacts = []
    for hres, hatom in hotspot_atoms:
        for bchain, bres, batom in binder_atoms:
            diff = hatom.coord - batom.coord
            dist = (diff @ diff) ** 0.5
            if dist <= cutoff:
                h_label = f"{target_chain_id}.{hres.id[1]}.{hatom.name.strip()}"
                b_label = f"{bchain.id}.{bres.id[1]}.{batom.name.strip()}"
                contacts.append((f"{h_label} - {b_label}", round(dist, 2)))

    return contacts


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Load config
    if not os.path.isfile(args.config):
        sys.exit(f"Error: config file '{args.config}' not found.")
    cfg = load_config(args.config)

    target_chain = cfg["chain"]
    hotspot_resnums = [int(r) for r in cfg["residues"].split(",")]
    cutoff = float(cfg["cutoff"])

    input_dir = args.input_dir or cfg.get("input_dir", "designs")
    output_dir = args.output_dir or cfg.get("output_dir", "designs_filtered")
    os.makedirs(output_dir, exist_ok=True)

    parser = PDBParser(QUIET=True)

    # Scan PDB files
    pdb_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".pdb"))
    if not pdb_files:
        sys.exit(f"No .pdb files found in '{input_dir}/'.")

    print(f"Scanning {len(pdb_files)} PDB files  "
          f"(chain={target_chain}, residues={hotspot_resnums}, cutoff={cutoff} A)")
    print()

    summary_rows = []
    passed = 0

    for filename in pdb_files:
        pdb_path = os.path.join(input_dir, filename)
        structure = parser.get_structure("model", pdb_path)
        contacts = find_contacts(structure, target_chain, hotspot_resnums, cutoff)

        if contacts:
            shutil.copy(pdb_path, os.path.join(output_dir, filename))
            passed += 1
            print(f"  PASS  {filename}  ({len(contacts)} contacts)")
            row = [filename, len(contacts)]
            for label, dist in contacts:
                row.extend([label, dist])
            summary_rows.append(row)
        else:
            print(f"  FAIL  {filename}")

    # Write CSV summary
    if summary_rows:
        max_contacts = max((len(row) - 2) // 2 for row in summary_rows)
        headers = ["filename", "contact_count"]
        for i in range(1, max_contacts + 1):
            headers.extend([f"contact_{i}", f"distance_{i}_A"])

        with open(args.summary, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(headers)
            writer.writerows(summary_rows)

    print()
    print(f"Done. {passed}/{len(pdb_files)} designs passed -> {output_dir}/")
    if summary_rows:
        print(f"Summary written to {args.summary}")


if __name__ == "__main__":
    main()
