"""
Step 1: Filter assembled ligands and compute Morgan fingerprints.

Input:  enamine_ligands.csv  (output of 0_build_enamine.py)
Output: <out_dir>/filtered_ligands.csv   — ligands passing drug-likeness filters
        <out_dir>/fingerprints.npy       — Morgan fp array (n_ligands × fp_size)

Usage:
    python3 1_fingerprint.py --input enamine_ligands.csv --out_dir fp_out/

    # With custom filters:
    python3 1_fingerprint.py --input enamine_ligands.csv --out_dir fp_out/ \\
        --mw_max 500 --logp_min 1 --logp_max 5 --hbd_max 3 --hba_max 8
"""
import argparse
import os
from time import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator


# ── Filters ──────────────────────────────────────────────────────────────────

def passes_filters(smi, mw_max, logp_min, logp_max, hbd_max, hba_max):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    # Exclude silicon-containing compounds (poor ADMET, common artefact)
    if any(a.GetAtomicNum() == 14 for a in mol.GetAtoms()):
        return False
    return (mw   <= mw_max and
            logp_min <= logp <= logp_max and
            hbd  <= hbd_max and
            hba  <= hba_max)


# ── Fingerprints ─────────────────────────────────────────────────────────────

def compute_fingerprints(smiles_list, fp_size=1024, radius=2):
    """Return (fp_array, valid_indices) for the given SMILES list."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    valid_idx = [i for i, m in enumerate(mols) if m is not None]
    valid_mols = [mols[i] for i in valid_idx]

    fp_array = np.zeros((len(valid_mols), fp_size), dtype=np.int8)
    t0 = time()
    for k, mol in enumerate(valid_mols):
        fp_array[k] = list(gen.GetFingerprint(mol))
        if k % 50_000 == 0 and k > 0:
            print(f"  fingerprints: {k:,}/{len(valid_mols):,}  ({time()-t0:.0f}s)")

    return fp_array, valid_idx


# ── Main ─────────────────────────────────────────────────────────────────────

def main(args):
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load ──
    print(f"\n── Loading ligands ──────────────────────────────────────────")
    df = pd.read_csv(args.input)
    print(f"  Total ligands: {len(df):,}")

    # ── Filter ──
    print(f"\n── Applying drug-likeness filters ───────────────────────────")
    print(f"  MW ≤ {args.mw_max}  |  logP {args.logp_min}–{args.logp_max}  "
          f"|  HBD ≤ {args.hbd_max}  |  HBA ≤ {args.hba_max}  |  no Si")
    mask = df['smiles'].apply(
        lambda s: passes_filters(s, args.mw_max, args.logp_min,
                                 args.logp_max, args.hbd_max, args.hba_max)
    )
    df_filtered = df[mask].reset_index(drop=True)
    print(f"  Passed: {len(df_filtered):,} / {len(df):,} "
          f"({100*len(df_filtered)/len(df):.1f}%)")

    out_csv = os.path.join(args.out_dir, 'filtered_ligands.csv')
    df_filtered.to_csv(out_csv, index=False)
    print(f"  Saved → {out_csv}")

    # ── Fingerprints ──
    print(f"\n── Computing Morgan fingerprints (r={args.radius}, size={args.fp_size}) ─")
    t0 = time()
    fp_array, valid_idx = compute_fingerprints(
        df_filtered['smiles'].tolist(), fp_size=args.fp_size, radius=args.radius
    )
    print(f"  Done: {fp_array.shape[0]:,} molecules in {time()-t0:.0f}s")

    # Drop rows whose SMILES failed to parse (should be near-zero after step 0)
    if len(valid_idx) < len(df_filtered):
        n_bad = len(df_filtered) - len(valid_idx)
        print(f"  [warn] {n_bad} molecules failed RDKit parse — dropped from fingerprint array")
        df_filtered = df_filtered.iloc[valid_idx].reset_index(drop=True)
        df_filtered.to_csv(out_csv, index=False)  # overwrite with clean set

    out_npy = os.path.join(args.out_dir, 'fingerprints.npy')
    np.save(out_npy, fp_array)
    print(f"  Saved → {out_npy}  shape={fp_array.shape}")

    print(f"\n── Done ─────────────────────────────────────────────────────")
    print(f"  Next step:")
    print(f"    python3 2_cluster_pick.py \\")
    print(f"      --fingerprints {out_npy} \\")
    print(f"      --ligands      {out_csv} \\")
    print(f"      --out_dir      <your_cluster_dir>/")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter assembled ligands and compute Morgan fingerprints.'
    )
    parser.add_argument('--input',    required=True,
                        help='Path to enamine_ligands.csv from step 0.')
    parser.add_argument('--out_dir',  default='fp_out',
                        help='Output directory (default: fp_out/).')
    # Filter params
    parser.add_argument('--mw_max',   type=float, default=450,
                        help='Maximum MW (default: 450).')
    parser.add_argument('--logp_min', type=float, default=2.0,
                        help='Minimum logP (default: 2.0).')
    parser.add_argument('--logp_max', type=float, default=4.0,
                        help='Maximum logP (default: 4.0).')
    parser.add_argument('--hbd_max',  type=int,   default=2,
                        help='Maximum H-bond donors (default: 2).')
    parser.add_argument('--hba_max',  type=int,   default=7,
                        help='Maximum H-bond acceptors (default: 7).')
    # Fingerprint params
    parser.add_argument('--fp_size',  type=int,   default=1024,
                        help='Morgan fingerprint bit size (default: 1024).')
    parser.add_argument('--radius',   type=int,   default=2,
                        help='Morgan radius (default: 2).')
    main(parser.parse_args())
