#!/usr/bin/env python3
"""
1b_reference_filter.py — Filter assembled enamine ligands by similarity to
                          known TgCPL lead compounds.

After 1_fingerprint.py has applied drug-likeness filters, this script further
narrows the library by Tanimoto similarity to the best TgCPL leads from the
paper (Cpd1, Cpd12, Cpd13, Cpd18).  Compounds with max Tanimoto >= --min_sim
(default 0.30) to any reference lead pass.

Optionally applies TgCPL-specific structural filters grounded in the SAR:
  --require_basic_n  : at least one basic sp3 N (engages TgCPL Glu75/Asp218)
  --require_scaffold : cyanotriazine-morpholine warhead must be present

Scientific rationale
--------------------
TgCPL's S2 pocket is narrow and predominantly acidic (Glu75, Asp218, Asp78),
while HsCPL's S2 is larger and hydrophobic.  High-selectivity leads carry a
basic amine at P2 that contacts these acidic residues.  Requiring a basic N
therefore enriches for TgCPL-selective chemotypes and strips compounds whose
P2 would fit HsCPL better.

Compounds are ranked by their max Tanimoto to any reference lead.  The
filtered fingerprint array is written alongside the filtered CSV so
2_cluster_pick.py can consume it directly without re-fingerprinting.

Usage
-----
  python 1b_reference_filter.py \\
      --input   fp_out/filtered_ligands.csv \\
      --fps     fp_out/fingerprints.npy \\
      --out_dir sim_out \\
      --min_sim 0.30 \\
      --require_basic_n \\
      --require_scaffold

Outputs (all in --out_dir)
--------------------------
  similar_ligands.csv          filtered ligand table (+ tanimoto_max column)
  similar_fps.npy              fingerprint rows aligned to similar_ligands.csv
  reference_filter_report.txt  stats per reference lead + filter tallies

Pipeline position
-----------------
  0_build_enamine.py      ->  enamine_ligands.csv
  1_fingerprint.py        ->  filtered_ligands.csv  fingerprints.npy
  1b_reference_filter.py  ->  similar_ligands.csv   similar_fps.npy    <- this
  2_cluster_pick.py       ->  picked_ligands.csv
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# ---------------------------------------------------------------------------
# Reference leads from: TgCPL Main Draft 6 / ShEPhERD-for-TgCPL leads.tsv
# Cpd18 is the best lead (207 nM TgCPL, 33x selectivity).
# Cpd12/13 have the best selectivity (cyclohexyl/cycloheptyl P3).
# Cpd1 is included for completeness (parent compound, lowest selectivity).
# ---------------------------------------------------------------------------
TGCPL_LEADS: dict[str, str] = {
    "Cpd1":  "CCN(CC)CCN(Cc1ccc(F)cc1)c1nc(C#N)nc(N2CCOCC2)n1",
    "Cpd12": "CCN(CC)CCN(c1nc(C#N)nc(N2CCOCC2)n1)C(c1ccc(F)cc1)C1CCCCC1",
    "Cpd13": "CCN(CC)CCN(c1nc(C#N)nc(N2CCOCC2)n1)C(c1ccc(F)cc1)C1CCCCCC1",
    "Cpd18": "N#Cc1nc(N2CCOCC2)nc(N(CCC2CCNCC2)C(c2ccc(F)cc2)C2CCCCC2)n1",
}

# SMARTS for structural filters
# Basic sp3 N: not aromatic, not amide, not sulfonamide, not guanidine
_BASIC_N_SMARTS       = "[NX3;H0,H1,H2;!$(NC=O);!$(NS=O);!$(NC=[NH]);!n]"
_CYANOTRIAZINE_SMARTS = "N#Cc1ncncn1"
_MORPHOLINE_SMARTS    = "N1CCOCC1"


def _morgan_fp(mol: Chem.Mol, radius: int, fp_size: int) -> DataStructs.ExplicitBitVect:
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_size)


def build_ref_fps(
    radius: int, fp_size: int
) -> dict[str, DataStructs.ExplicitBitVect]:
    fps: dict[str, DataStructs.ExplicitBitVect] = {}
    for name, smi in TGCPL_LEADS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  WARNING: could not parse reference SMILES for {name}", file=sys.stderr)
            continue
        fps[name] = _morgan_fp(mol, radius, fp_size)
    return fps


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="filtered_ligands.csv from 1_fingerprint.py")
    p.add_argument("--fps",     required=True,
                   help="fingerprints.npy from 1_fingerprint.py (rows aligned to --input)")
    p.add_argument("--out_dir", default="sim_out",
                   help="output directory (default: sim_out)")
    p.add_argument("--min_sim", type=float, default=0.30,
                   help="minimum Tanimoto to any reference lead (default: 0.30)")
    p.add_argument("--radius",  type=int,   default=2,
                   help="Morgan FP radius (default: 2; must match 1_fingerprint.py)")
    p.add_argument("--fp_size", type=int,   default=1024,
                   help="FP bit size (default: 1024; must match 1_fingerprint.py)")
    p.add_argument("--require_basic_n",  action="store_true",
                   help="require at least one basic sp3 N (enriches for TgCPL S2 contacts)")
    p.add_argument("--require_scaffold", action="store_true",
                   help="require the cyanotriazine-morpholine warhead (N#Cc1ncncn1 + morpholine)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()

    # -----------------------------------------------------------------------
    # Load inputs
    # -----------------------------------------------------------------------
    print(f"Loading {args.input} ...")
    df = pd.read_csv(args.input)
    print(f"  {len(df):,} drug-like compounds")

    print(f"Loading fingerprints from {args.fps} ...")
    fps_arr = np.load(args.fps)
    if fps_arr.shape[0] != len(df):
        sys.exit(
            f"ERROR: fingerprint rows ({fps_arr.shape[0]}) != CSV rows ({len(df)}). "
            "Re-run 1_fingerprint.py to regenerate aligned outputs."
        )

    # -----------------------------------------------------------------------
    # Build reference fingerprints
    # -----------------------------------------------------------------------
    print("Building reference fingerprints ...")
    ref_fps = build_ref_fps(args.radius, args.fp_size)
    if not ref_fps:
        sys.exit("ERROR: no valid reference fingerprints.")
    print(f"  References: {', '.join(ref_fps)}")
    ref_bv_list = list(ref_fps.values())
    ref_names   = list(ref_fps.keys())

    # -----------------------------------------------------------------------
    # Compile substructure patterns for optional filters
    # -----------------------------------------------------------------------
    basic_n_patt = Chem.MolFromSmarts(_BASIC_N_SMARTS)
    cn_patt      = Chem.MolFromSmarts(_CYANOTRIAZINE_SMARTS)
    morph_patt   = Chem.MolFromSmarts(_MORPHOLINE_SMARTS)

    # -----------------------------------------------------------------------
    # Score every compound
    # -----------------------------------------------------------------------
    print(f"Scoring {len(df):,} compounds (min_sim={args.min_sim}) ...")

    tanimoto_max     = np.full(len(df), -1.0, dtype=np.float32)
    tanimoto_per_ref = {n: np.full(len(df), -1.0, dtype=np.float32) for n in ref_names}
    basic_n_flag     = np.zeros(len(df), dtype=bool)
    scaffold_flag    = np.zeros(len(df), dtype=bool)

    n_invalid = 0
    smiles_col = "smiles" if "smiles" in df.columns else df.columns[1]

    for i, smi in enumerate(df[smiles_col]):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            n_invalid += 1
            continue

        fp   = _morgan_fp(mol, args.radius, args.fp_size)
        sims = DataStructs.BulkTanimotoSimilarity(fp, ref_bv_list)

        for j, name in enumerate(ref_names):
            tanimoto_per_ref[name][i] = sims[j]
        tanimoto_max[i] = max(sims)

        if args.require_basic_n:
            basic_n_flag[i] = mol.HasSubstructMatch(basic_n_patt)
        if args.require_scaffold:
            scaffold_flag[i] = (mol.HasSubstructMatch(cn_patt) and
                                mol.HasSubstructMatch(morph_patt))

        if (i + 1) % 100_000 == 0:
            elapsed = time.perf_counter() - t0
            rate    = (i + 1) / elapsed
            remaining = (len(df) - i - 1) / rate
            print(f"  {i+1:>8,}/{len(df):,}  "
                  f"elapsed {elapsed:.0f}s  "
                  f"remaining ~{remaining:.0f}s  "
                  f"({rate:.0f} cpd/s)", flush=True)

    # -----------------------------------------------------------------------
    # Apply filters
    # -----------------------------------------------------------------------
    keep = tanimoto_max >= args.min_sim
    n_after_sim = keep.sum()

    n_after_basic_n  = n_after_sim
    n_after_scaffold = n_after_sim

    if args.require_basic_n:
        keep &= basic_n_flag
        n_after_basic_n = keep.sum()

    if args.require_scaffold:
        keep &= scaffold_flag
        n_after_scaffold = keep.sum()

    # -----------------------------------------------------------------------
    # Build output dataframe
    # -----------------------------------------------------------------------
    df_out = df[keep].copy()
    df_out["tanimoto_max"] = tanimoto_max[keep]
    for name in ref_names:
        df_out[f"tanimoto_{name}"] = tanimoto_per_ref[name][keep]
    df_out.sort_values("tanimoto_max", ascending=False, inplace=True)
    df_out.reset_index(drop=True, inplace=True)

    fps_out = fps_arr[keep]

    # -----------------------------------------------------------------------
    # Write outputs
    # -----------------------------------------------------------------------
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "similar_ligands.csv"
    npy_path = out_dir / "similar_fps.npy"
    rpt_path = out_dir / "reference_filter_report.txt"

    df_out.to_csv(csv_path, index=False)
    np.save(npy_path, fps_out)

    elapsed = time.perf_counter() - t0

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    report_lines = [
        "=" * 60,
        "TgCPL reference similarity filter report",
        "=" * 60,
        f"Input compounds (drug-like):  {len(df):>10,}",
        f"Invalid SMILES (skipped):     {n_invalid:>10,}",
        "",
        "Filter funnel:",
        f"  Tanimoto >= {args.min_sim:.2f} to any ref: {n_after_sim:>10,}",
    ]
    if args.require_basic_n:
        report_lines.append(
            f"  + basic N required:          {n_after_basic_n:>10,}"
        )
    if args.require_scaffold:
        report_lines.append(
            f"  + scaffold required:         {n_after_scaffold:>10,}"
        )
    report_lines += [
        f"  PASSED (written to CSV):     {len(df_out):>10,}",
        f"  Retention:                   {100*len(df_out)/max(len(df),1):.1f}%",
        "",
        "Per-reference similarity stats (over passing compounds):",
    ]
    for name in ref_names:
        col = df_out[f"tanimoto_{name}"]
        report_lines.append(
            f"  {name:6s}  mean={col.mean():.3f}  "
            f"median={col.median():.3f}  max={col.max():.3f}"
        )
    report_lines += [
        "",
        "Top 10 by tanimoto_max:",
    ]
    top10_cols = ["smiles", "tanimoto_max"] + [f"tanimoto_{n}" for n in ref_names]
    top10_cols = [c for c in top10_cols if c in df_out.columns]
    report_lines.append(df_out.head(10)[top10_cols].to_string(index=False))
    report_lines += [
        "",
        f"Elapsed: {elapsed:.1f}s",
        "",
        "Output files:",
        f"  {csv_path}",
        f"  {npy_path}",
        "=" * 60,
    ]

    report_text = "\n".join(report_lines)
    print()
    print(report_text)
    rpt_path.write_text(report_text)

    print(f"\nDone. {len(df_out):,} compounds written to {out_dir}/")
    print(f"Pass these to 2_cluster_pick.py:")
    print(f"  python 2_cluster_pick.py \\")
    print(f"      --fingerprints {npy_path} \\")
    print(f"      --ligands      {csv_path} \\")
    print(f"      --out_dir      cluster_out")


if __name__ == "__main__":
    main()
