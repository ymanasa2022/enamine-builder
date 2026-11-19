"""
build_enamine.py  —  Combinatorial triazine library builder (optimised)

Scaffold:  2-CN-4,6-dichloro-1,3,5-triazine   N#CC1=NC(Cl)=NC(Cl)=N1
           P1  = variable (secondary amine OR alcohol) via SNAr at first Cl
           P2/P3 = secondary amine arm  (R'–CH₂–NH–CH₂–R)  via SNAr at second Cl

Building-block pools (from ENAMINE SDF files):
  --alcohols        → alcohols / OH sources
  --aldehydes       → aldehyde sources  (required)
  --amines          → primary + secondary amines  (required)
  --ester_carboxyl  → combined ester+carboxyl file (optional)
  --halides         → alkyl halides for SN2 (optional)

Speed notes
  - All SMARTS reactions pre-compiled at module level (not per-call)
  - Molecules parsed once per pool into (key, mol) pairs
  - Cross-product loops parallelised with multiprocessing (--workers N)
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFingerprintGenerator
import argparse
import pandas as pd
import os
import datetime
import glob as _glob
from multiprocessing import Pool, cpu_count


# ---------------------------------------------------------------------------
# TgCPL reference P2/P3 secondary amines (extracted from best leads).
# Used by --prefilter_p23_sim to rank the P2/P3 pool before assembly.
# Cpd18 is the best lead (207 nM TgCPL, 33x sel) — piperidinylethyl P2 +
# cyclohexyl-4-fluorobenzyl P3.  Cpd12/13 have optimal selectivity.
# The direct secondary amine that attaches to triazine via SNAr is used.
# ---------------------------------------------------------------------------
_P23_REFERENCE_SMILES = [
    # Cpd18 P2/P3 amine: piperidinylethyl × cyclohexyl-4-fluorobenzyl
    "[NH](CCC1CCNCC1)C(c1ccc(F)cc1)C1CCCCC1",
    # Cpd12 P2/P3 amine: diethylaminoethyl × cyclohexyl-4-fluorobenzyl
    "[NH](CCN(CC)CC)C(c1ccc(F)cc1)C1CCCCC1",
    # Cpd13 P2/P3 amine: diethylaminoethyl × cycloheptyl-4-fluorobenzyl
    "[NH](CCN(CC)CC)C(c1ccc(F)cc1)C1CCCCCC1",
]


def _prefilter_p23_by_similarity(p23_dict, min_sim, radius=2, fp_size=1024):
    """Return a subset of p23_dict keeping entries with Tanimoto >= min_sim
    to any of the _P23_REFERENCE_SMILES.

    Rational basis: TgCPL's S2 pocket is acidic (Glu75/Asp218/Asp78).
    High-selectivity leads carry a basic amine P2 that contacts these residues,
    plus a bulky cycloalkyl P3 that sterically forces selectivity over HsCPL's
    wider, hydrophobic S2.  Filtering the 10+ M P2/P3 pool to fragments that
    resemble these patterns before assembly reduces the combinatorial explosion
    by 10-50x with minimal loss of on-target chemistry.

    Args:
        p23_dict : {name: smiles} P2/P3 secondary amine pool
        min_sim  : Tanimoto threshold (0.15-0.25 is a practical range)
        radius   : Morgan FP radius (must match downstream fingerprinting)
        fp_size  : Morgan FP bit width

    Returns:
        kept : {name: smiles} filtered subset
    """
    import time as _time
    t0 = _time.perf_counter()

    # Build reference fps
    ref_fps = []
    for smi in _P23_REFERENCE_SMILES:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            ref_fps.append(rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size).GetFingerprint(mol))
    if not ref_fps:
        print("  WARNING: could not parse any reference SMILES; skipping prefilter.")
        return p23_dict

    kept = {}
    n = len(p23_dict)
    for i, (key, smi) in enumerate(p23_dict.items()):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp   = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size).GetFingerprint(mol)
        sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        if max(sims) >= min_sim:
            kept[key] = smi
        if (i + 1) % 500_000 == 0:
            elapsed = _time.perf_counter() - t0
            print(f"    prefilter: {i+1:>8,}/{n:,}  kept so far: {len(kept):,}  "
                  f"({elapsed:.0f}s)", flush=True)

    elapsed = _time.perf_counter() - t0
    print(f"    prefilter done in {elapsed:.1f}s: "
          f"{n:,} → {len(kept):,} P2/P3 entries "
          f"(sim >= {min_sim}, {100*len(kept)/max(n,1):.1f}% retained)")
    return kept

now = datetime.datetime.now()


# ── Pre-compiled reactions (module level — pay the cost once) ─────────────────

_RXN_ESTER_TO_ALD    = AllChem.ReactionFromSmarts('[C:1](=[O])[O][C]>>[C:1]O')
_RXN_CARBOXYL_TO_ALD = AllChem.ReactionFromSmarts('[C:0][C:1](=[O])[OH]>>[C:0][CH:1]=O')
_RXN_OH_TO_ALD       = AllChem.ReactionFromSmarts('[C:0][CH2:1][OH:2]>>[C:0][CH:1]=[O:2]')
_RXN_ALD_TO_OH       = AllChem.ReactionFromSmarts('[C:0][CH:1]=[O:2]>>[C:0][CH2:1][OH:2]')
_RXN_ESTER_TO_OH     = AllChem.ReactionFromSmarts('[C:0][C:1](=[O])[O][C]>>[C:0][CH2:1][OH]')
_RXN_CARBOXYL_TO_OH  = AllChem.ReactionFromSmarts('[C:0][C:1](=[O])[OH]>>[C:0][CH2:1][OH]')
_RXN_NITRILE_TO_NH2  = AllChem.ReactionFromSmarts('[C:1]#[N:2]>>[C:1][NH2:2]')
_RXN_IMINE_TO_NH2    = AllChem.ReactionFromSmarts('[C:1]=[NH:2]>>[C:1][NH2:2]')
_RXN_REDUCTIVE_AM    = AllChem.ReactionFromSmarts('[NH2:1].[CH:2]=[O]>>[NH:1][CH2:2]')
_RXN_SN2             = AllChem.ReactionFromSmarts('[NH2:1].[C:2][Cl,Br,I]>>[NH:1][C:2]')
_RXN_SNAR_N          = AllChem.ReactionFromSmarts('[c:1][Cl].[NH:3]>>[c:1][N:3]')
_RXN_SNAR_O          = AllChem.ReactionFromSmarts('[c:1][Cl].[C:2][OH]>>[c:1][O][C:2]')

_PAT_CL_AROM    = Chem.MolFromSmarts('[c][Cl]')
_PAT_PRIMARY_OH = Chem.MolFromSmarts('[C][CH2][OH]')
_PAT_PRIM_AMINE = Chem.MolFromSmarts('[C][NH2]')
_PAT_SEC_AMINE  = Chem.MolFromSmarts('[C][NH][C]')
_PAT_NITRILE    = Chem.MolFromSmarts('[C]#[N]')
_PAT_IMINE      = Chem.MolFromSmarts('[C]=[NH]')
_PAT_HALIDE     = Chem.MolFromSmarts('[C][F,Cl,Br,I]')

CORE_DICHLORO = 'N#CC1=NC(Cl)=NC(Cl)=N1'


# ── SDF loader ───────────────────────────────────────────────────────────────

def _mol_id(mol, idx):
    for prop in ('CatalogID', 'ID', 'Code'):
        if mol.HasProp(prop):
            v = mol.GetProp(prop).strip()
            if v:
                return v
    name = mol.GetProp('_Name').strip() if mol.HasProp('_Name') else ''
    return name if name else f'mol_{idx}'


def load_sdf(path):
    suppl = Chem.SDMolSupplier(path, sanitize=True, removeHs=True)
    records = []
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        records.append({'synonyms': _mol_id(mol, i), 'smiles': Chem.MolToSmiles(mol)})
    return pd.DataFrame(records)


def load_block(enamine_dir, prefix, label):
    paths = sorted(_glob.glob(os.path.join(enamine_dir, f'{prefix}*.sdf')))
    if not paths:
        raise FileNotFoundError(f"No SDF files for {label} (prefix '{prefix}' in {enamine_dir})")
    frames = []
    for p in paths:
        df = load_sdf(p)
        print(f"  {os.path.basename(p):35s}  {len(df):>7,} molecules")
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    print(f"  → {label} total: {len(combined):,}")
    return combined


def _load_optional(enamine_dir, prefix, label):
    if not prefix:
        print(f"  [skip] {label} — not provided")
        return pd.DataFrame(columns=['synonyms', 'smiles'])
    try:
        return load_block(enamine_dir, prefix, label)
    except FileNotFoundError as e:
        print(f"  [warn] {e}")
        return pd.DataFrame(columns=['synonyms', 'smiles'])


# ── Property helpers ─────────────────────────────────────────────────────────

def compute_props(smi, name=''):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return {
        'name':      name,
        'smiles':    smi,
        'MW':        round(Descriptors.MolWt(mol), 2),
        'logP':      round(Descriptors.MolLogP(mol), 2),
        'HBD':       Descriptors.NumHDonors(mol),
        'HBA':       Descriptors.NumHAcceptors(mol),
        'RotBonds':  rdMolDescriptors.CalcNumRotatableBonds(mol),
        'TPSA':      round(rdMolDescriptors.CalcTPSA(mol), 1),
        'RingCount': mol.GetRingInfo().NumRings(),
    }


def log_pool(pool_dict, label, out_dir):
    records = [compute_props(smi, name) for name, smi in pool_dict.items() if smi]
    records = [r for r in records if r]
    if not records:
        print(f"  [warn] {label}: nothing to log")
        return
    df = pd.DataFrame(records)
    path = os.path.join(out_dir, f'props_{label}.csv')
    df.to_csv(path, index=False)
    print(f"  {label}: {len(df):,} entries → {path}")
    print(f"    MW {df['MW'].min():.0f}–{df['MW'].max():.0f}  "
          f"logP {df['logP'].min():.1f}–{df['logP'].max():.1f}  "
          f"HBD {df['HBD'].min()}–{df['HBD'].max()}  "
          f"RotBonds {df['RotBonds'].min()}–{df['RotBonds'].max()}")


# ── Block-level filters ──────────────────────────────────────────────────────

def filter_p3_block(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw  = Descriptors.MolWt(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rng = mol.GetRingInfo().NumRings()
    return (mw > 100 or rng >= 1) and mw <= 200 and rot <= 3


def filter_p2_block(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw  = Descriptors.MolWt(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return (hbd + hba) >= 1 and mw <= 200 and rot <= 3


def filter_p1_block(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw  = Descriptors.MolWt(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return mw <= 150 and rot <= 3


def filter_post_p23(smi):
    """Post-synthesis filter on assembled P2/P3 secondary amines.
    Keeps products that are small enough to leave room for the triazine
    scaffold (~150 Da) and P1 block (~100 Da) in the final molecule.
    Target final MW budget: ~500 Da → P2/P3 arm ≤ 300 Da.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw  = Descriptors.MolWt(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return mw <= 300 and rot <= 5


def print_filter_summary(label, before, after):
    pct = 100 * after / before if before > 0 else 0
    print(f"  {label}: {before:,} → {after:,} ({pct:.1f}% kept)")


# ── Mol-level conversions (accept mol object, return SMILES or None) ──────────

def _apply_rxn(rxn, mol):
    """Run a single-reactant reaction; return first valid product SMILES or None."""
    try:
        prods = rxn.RunReactants((mol,))
        if prods:
            p = Chem.MolFromSmiles(Chem.MolToSmiles(prods[0][0]))
            if p:
                return Chem.MolToSmiles(p)
    except Exception:
        pass
    return None


def get_aldehyde(mol):
    # Already an aldehyde? (–CHO)
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[CH]=O')):
        return Chem.MolToSmiles(mol)
    # Ester → intermediate OH, then oxidise
    oh_smi = _apply_rxn(_RXN_ESTER_TO_ALD, mol)
    if oh_smi:
        oh_mol = Chem.MolFromSmiles(oh_smi)
        result = _apply_rxn(_RXN_OH_TO_ALD, oh_mol) if oh_mol else None
        if result:
            return result
    # Carboxyl → aldehyde directly
    result = _apply_rxn(_RXN_CARBOXYL_TO_ALD, mol)
    if result:
        return result
    # Primary alcohol → aldehyde
    return _apply_rxn(_RXN_OH_TO_ALD, mol)


def get_alcohol(mol):
    if mol.HasSubstructMatch(_PAT_PRIMARY_OH):
        return Chem.MolToSmiles(mol)
    result = _apply_rxn(_RXN_ALD_TO_OH, mol)
    if result:
        return result
    result = _apply_rxn(_RXN_ESTER_TO_OH, mol)
    if result:
        return result
    return _apply_rxn(_RXN_CARBOXYL_TO_OH, mol)


def classify_and_process_amine(smi):
    """Return (processed_smiles, 'primary'|'secondary'|'other')."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'invalid'
    # Convert nitrile / imine precursors to primary amine
    if mol.HasSubstructMatch(_PAT_NITRILE):
        new_smi = _apply_rxn(_RXN_NITRILE_TO_NH2, mol)
        if new_smi:
            mol = Chem.MolFromSmiles(new_smi)
    if mol.HasSubstructMatch(_PAT_IMINE):
        new_smi = _apply_rxn(_RXN_IMINE_TO_NH2, mol)
        if new_smi:
            mol = Chem.MolFromSmiles(new_smi)
    out = Chem.MolToSmiles(mol)
    if mol.HasSubstructMatch(_PAT_PRIM_AMINE):
        return out, 'primary'
    if mol.HasSubstructMatch(_PAT_SEC_AMINE):
        return out, 'secondary'
    return out, 'other'


# ── Parallel worker functions (must be top-level for pickling) ────────────────
#
# Design: shared lists (aldehydes, halides, p23 amines) are sent to each worker
# process ONCE via the pool initializer, stored in _G_SHARED.  Individual tasks
# are just (key, smiles) pairs — tiny — so IPC overhead is minimal even for
# millions of cross-product pairs.

_G_SHARED: dict = {}


def _init_worker(shared: dict):
    """Pool initializer: load shared data into the worker-process global."""
    global _G_SHARED
    _G_SHARED = shared


def _worker_reductive_am(amine_kv):
    """One primary amine × all aldehydes → [(name, smiles), …]
    Applies filter_post_p23 before returning to avoid sending large
    volumes of data back over the multiprocessing pipe.
    """
    amine_key, amine_smi = amine_kv
    amine = Chem.MolFromSmiles(amine_smi)
    if amine is None:
        return []
    out = []
    for ald_key, ald_smi in _G_SHARED['ald_items']:
        ald = Chem.MolFromSmiles(ald_smi)
        if ald is None:
            continue
        try:
            prods = _RXN_REDUCTIVE_AM.RunReactants((amine, ald))
            if prods:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(prods[0][0]))
                if mol:
                    smi = Chem.MolToSmiles(mol)
                    if filter_post_p23(smi):
                        out.append((f'{amine_key}_{ald_key}', smi))
        except Exception:
            pass
    return out


def _worker_sn2(amine_kv):
    """One primary amine × all halides → [(name, smiles), …]
    Applies filter_post_p23 before returning to avoid sending large
    volumes of data back over the multiprocessing pipe.
    """
    amine_key, amine_smi = amine_kv
    amine = Chem.MolFromSmiles(amine_smi)
    if amine is None:
        return []
    out = []
    for hal_key, hal_smi in _G_SHARED['hal_items']:
        hal = Chem.MolFromSmiles(hal_smi)
        if hal is None:
            continue
        try:
            prods = _RXN_SN2.RunReactants((amine, hal))
            if prods:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(prods[0][0]))
                if mol:
                    smi = Chem.MolToSmiles(mol)
                    if filter_post_p23(smi):
                        out.append((f'{amine_key}_{hal_key}', smi))
        except Exception:
            pass
    return out


def _worker_assemble_amine_p1(p1_kv):
    """One P1 secondary amine × all P23 amines → [row_dict, …]"""
    p1_key, p1_smi = p1_kv
    p1_mol = Chem.MolFromSmiles(p1_smi)
    if p1_mol is None:
        return []
    mono = None
    try:
        core = Chem.MolFromSmiles(CORE_DICHLORO)
        prods = _RXN_SNAR_N.RunReactants((core, p1_mol))
        for ps in prods:
            m = Chem.MolFromSmiles(Chem.MolToSmiles(ps[0]))
            if m and m.HasSubstructMatch(_PAT_CL_AROM):
                mono = Chem.MolToSmiles(m)
                break
    except Exception:
        pass
    if mono is None:
        return []
    mono_mol = Chem.MolFromSmiles(mono)
    out = []
    for p23_key, p23_smi in _G_SHARED['p23_items']:
        p23_mol = Chem.MolFromSmiles(p23_smi)
        if p23_mol is None:
            continue
        try:
            prods = _RXN_SNAR_N.RunReactants((mono_mol, p23_mol))
            if prods:
                full = Chem.MolFromSmiles(Chem.MolToSmiles(prods[0][0]))
                if full:
                    out.append({
                        'synonyms':   f'{p1_key}::{p23_key}',
                        'p1_name':    p1_key,
                        'p1_smiles':  p1_smi,
                        'p1_type':    'amine',
                        'p23_name':   p23_key,
                        'p23_smiles': p23_smi,
                        'smiles':     Chem.MolToSmiles(full),
                    })
        except Exception:
            pass
    return out


def _worker_assemble_alcohol_p1(p1_kv):
    """One P1 alcohol × all P23 amines → [row_dict, …]"""
    p1_key, p1_smi = p1_kv
    p1_mol = Chem.MolFromSmiles(p1_smi)
    if p1_mol is None:
        return []
    mono = None
    try:
        core = Chem.MolFromSmiles(CORE_DICHLORO)
        prods = _RXN_SNAR_O.RunReactants((core, p1_mol))
        for ps in prods:
            m = Chem.MolFromSmiles(Chem.MolToSmiles(ps[0]))
            if m and m.HasSubstructMatch(_PAT_CL_AROM):
                mono = Chem.MolToSmiles(m)
                break
    except Exception:
        pass
    if mono is None:
        return []
    mono_mol = Chem.MolFromSmiles(mono)
    out = []
    for p23_key, p23_smi in _G_SHARED['p23_items']:
        p23_mol = Chem.MolFromSmiles(p23_smi)
        if p23_mol is None:
            continue
        try:
            prods = _RXN_SNAR_N.RunReactants((mono_mol, p23_mol))
            if prods:
                full = Chem.MolFromSmiles(Chem.MolToSmiles(prods[0][0]))
                if full:
                    out.append({
                        'synonyms':   f'{p1_key}::{p23_key}',
                        'p1_name':    p1_key,
                        'p1_smiles':  p1_smi,
                        'p1_type':    'alcohol',
                        'p23_name':   p23_key,
                        'p23_smiles': p23_smi,
                        'smiles':     Chem.MolToSmiles(full),
                    })
        except Exception:
            pass
    return out


def _parallel(worker_fn, task_list, n_workers, chunksize=20, shared=None):
    """Run worker_fn over task_list using a process pool; flatten results.

    shared: dict passed to each worker process once via the pool initializer.
            Use this for large shared lists (aldehydes, p23 amines, etc.) to
            avoid re-pickling them for every task.
    """
    if n_workers == 1:
        if shared:
            _init_worker(shared)
        results = [worker_fn(t) for t in task_list]
    else:
        init_args = (shared,) if shared else ({},)
        with Pool(n_workers, initializer=_init_worker, initargs=init_args) as pool:
            results = pool.map(worker_fn, task_list, chunksize=chunksize)
    return [item for batch in results for item in batch]


# ── Main ─────────────────────────────────────────────────────────────────────

def main(args):
    out_dir  = args.out_dir
    workers  = args.workers
    os.makedirs(out_dir, exist_ok=True)
    print(f"\nUsing {workers} worker(s)")

    # ── Resume path: load pre-built pools and skip to assembly ───────────────
    if args.skip_synthesis:
        def _load_pool(label):
            path = os.path.join(out_dir, f'props_{label}.csv')
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"--skip_synthesis requires {path} from a previous run. "
                    f"Run without --skip_synthesis first.")
            df = pd.read_csv(path)
            return dict(zip(df['name'], df['smiles']))

        print("\n── Skipping synthesis — loading saved pools from", out_dir, "────────")
        sec_amines_p1   = _load_pool('secondary_amines_p1')
        alcohols        = _load_pool('alcohols')
        p23_sec_amines  = _load_pool('p23_secondary_amines')
        print(f"  sec_amines_p1 : {len(sec_amines_p1):,}")
        print(f"  alcohols      : {len(alcohols):,}")
        print(f"  p23_sec_amines: {len(p23_sec_amines):,}")

        if args.prefilter_p23_sim is not None:
            print(f"\n── P2/P3 reference pre-filter (sim >= {args.prefilter_p23_sim}) ────────────")
            p23_sec_amines = _prefilter_p23_by_similarity(
                p23_sec_amines, args.prefilter_p23_sim
            )

        # Jump straight to assembly (step 7)
        print("\n── Assembly ─────────────────────────────────────────────────────────")
        p23_items = list(p23_sec_amines.items())
        result_rows = []
        assemble_shared = {'p23_items': p23_items}

        a_tasks = list(sec_amines_p1.items())
        print(f"  Pipeline A (P1=amine): {len(sec_amines_p1):,} × {len(p23_sec_amines):,} "
              f"= {len(sec_amines_p1)*len(p23_sec_amines):,} pairs …")
        result_rows += _parallel(_worker_assemble_amine_p1, a_tasks, workers,
                                 shared=assemble_shared)
        print(f"    → {len(result_rows):,} products")

        if getattr(args, 'skip_alcohol_p1', False):
            print("  Pipeline B (P1=alcohol): SKIPPED (--skip_alcohol_p1)")
        else:
            b_tasks  = list(alcohols.items())
            before_b = len(result_rows)
            print(f"  Pipeline B (P1=alcohol): {len(alcohols):,} × {len(p23_sec_amines):,} "
                  f"= {len(alcohols)*len(p23_sec_amines):,} pairs …")
            result_rows += _parallel(_worker_assemble_alcohol_p1, b_tasks, workers,
                                     shared=assemble_shared)
            print(f"    → {len(result_rows)-before_b:,} products")

        result_df = pd.DataFrame(result_rows)
        print(f"\n  Total assembled: {len(result_df):,}")

        out_path = os.path.join(args.enamine_dir, 'enamine_ligands.csv')
        if os.path.exists(out_path):
            out_path = os.path.join(args.enamine_dir,
                                    f'enamine_ligands_{now.strftime("%Y%m%d_%H%M%S")}.csv')
        result_df[['synonyms', 'smiles']].to_csv(out_path, index=False)
        print(f"\nSaved (minimal) → {out_path}")
        full_path = out_path.replace('.csv', '_full.csv')
        result_df.to_csv(full_path, index=False)
        print(f"Saved (full)    → {full_path}")
        return

    # ── 1. Load SDF files ────────────────────────────────────────────────────
    print("\n── Loading SDF files ───────────────────────────────────────────────")
    df_oh    = load_block(args.enamine_dir, args.alcohols,  'alcohols')
    df_ald   = load_block(args.enamine_dir, args.aldehydes, 'aldehydes')
    df_amine = load_block(args.enamine_dir, args.amines,    'amines')

    df_ec     = _load_optional(args.enamine_dir, args.ester_carboxyl, 'esters+carboxyls')
    df_halide = _load_optional(args.enamine_dir, args.halides,         'halides')

    # ── 2. Build aldehyde pool ───────────────────────────────────────────────
    print("\n── Aldehyde pool ───────────────────────────────────────────────────")
    ald_src = pd.concat([df_ald, df_oh, df_ec], ignore_index=True)
    ald_src = ald_src[ald_src['smiles'].apply(filter_p3_block)].reset_index(drop=True)
    print(f"  Sources after P3 filter: {len(ald_src):,}")

    aldehydes = {}
    for _, row in ald_src.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol:
            result = get_aldehyde(mol)
            if result:
                aldehydes[row['synonyms']] = result
    print(f"  Converted to aldehydes: {len(aldehydes):,}")
    log_pool(aldehydes, 'aldehydes', out_dir)

    # ── 3. Build alcohol pool (P1 O-ether) ───────────────────────────────────
    print("\n── Alcohol pool ────────────────────────────────────────────────────")
    oh_src = pd.concat([df_oh, df_ald, df_ec], ignore_index=True)
    oh_src = oh_src[oh_src['smiles'].apply(filter_p1_block)].reset_index(drop=True)
    print(f"  Sources after P1 filter: {len(oh_src):,}")

    alcohols = {}
    for _, row in oh_src.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol:
            result = get_alcohol(mol)
            if result:
                alcohols[row['synonyms']] = result
    print(f"  Converted to alcohols: {len(alcohols):,}")
    log_pool(alcohols, 'alcohols', out_dir)

    # ── 4. Build halide pool (optional) ─────────────────────────────────────
    print("\n── Halide pool ─────────────────────────────────────────────────────")
    halides = {}
    if not df_halide.empty:
        for _, row in df_halide[df_halide['smiles'].apply(filter_p3_block)].iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol and mol.HasSubstructMatch(_PAT_HALIDE):
                halides[row['synonyms']] = Chem.MolToSmiles(mol)
        print(f"  Valid halides: {len(halides):,}")
        log_pool(halides, 'halides', out_dir)
    else:
        print("  [skip] no halide SDF — SN2 path skipped")

    # ── 5. Classify amines ───────────────────────────────────────────────────
    print("\n── Amine classification ────────────────────────────────────────────")
    primary_amines   = {}
    secondary_amines = {}
    for _, row in df_amine.iterrows():
        smi, aclass = classify_and_process_amine(row['smiles'])
        if aclass == 'primary':
            primary_amines[row['synonyms']] = smi
        elif aclass == 'secondary':
            secondary_amines[row['synonyms']] = smi

    primary_amines   = {k: v for k, v in primary_amines.items()   if filter_p2_block(v)}
    sec_amines_p1    = {k: v for k, v in secondary_amines.items() if filter_p1_block(v)}
    sec_amines_p23   = {k: v for k, v in secondary_amines.items() if filter_p2_block(v)}

    print(f"  Primary amines (P2/P3-ready): {len(primary_amines):,}")
    print(f"  Secondary amines for P1:      {len(sec_amines_p1):,}")
    print(f"  Secondary amines for P2/P3:   {len(sec_amines_p23):,}")
    log_pool(primary_amines, 'primary_amines',      out_dir)
    log_pool(sec_amines_p1,  'secondary_amines_p1', out_dir)

    # ── 6. Synthesise P2/P3 secondary amines (parallelised) ──────────────────
    print("\n── P2/P3 secondary amine synthesis ─────────────────────────────────")
    p23_sec_amines = {}

    ald_items = list(aldehydes.items())
    hal_items = list(halides.items())

    # Reductive amination  — shared list sent once per worker via initializer
    ra_tasks = list(primary_amines.items())
    print(f"  Reductive amination: {len(primary_amines):,} × {len(aldehydes):,} "
          f"= {len(primary_amines)*len(aldehydes):,} pairs …")
    for name, smi in _parallel(_worker_reductive_am, ra_tasks, workers,
                                shared={'ald_items': ald_items}):
        p23_sec_amines[name] = smi
    print(f"    → {len(p23_sec_amines):,} products")

    # SN2
    if halides:
        sn2_tasks = list(primary_amines.items())
        before = len(p23_sec_amines)
        print(f"  SN2: {len(primary_amines):,} × {len(halides):,} "
              f"= {len(primary_amines)*len(halides):,} pairs …")
        for name, smi in _parallel(_worker_sn2, sn2_tasks, workers,
                                    shared={'hal_items': hal_items}):
            p23_sec_amines[name] = smi
        print(f"    → {len(p23_sec_amines)-before:,} new products")

    # Direct secondary amines
    for k, v in sec_amines_p23.items():
        p23_sec_amines[f'direct_{k}'] = v
    print(f"  Total P2/P3 pool (incl. direct): {len(p23_sec_amines):,}")

    # Post-synthesis filter: MW ≤ 300, RotBonds ≤ 5
    before_post = len(p23_sec_amines)
    p23_sec_amines = {k: v for k, v in p23_sec_amines.items() if filter_post_p23(v)}
    print_filter_summary("  Post-synthesis P2/P3 filter (MW≤300, RotBonds≤5)",
                         before_post, len(p23_sec_amines))
    log_pool(p23_sec_amines, 'p23_secondary_amines', out_dir)

    # ── 6b. Optional P2/P3 similarity pre-filter ─────────────────────────────
    if args.prefilter_p23_sim is not None:
        print(f"\n── P2/P3 reference pre-filter (sim >= {args.prefilter_p23_sim}) ────────────")
        p23_sec_amines = _prefilter_p23_by_similarity(
            p23_sec_amines, args.prefilter_p23_sim
        )

    # ── 7. Assembly (parallelised) ────────────────────────────────────────────
    print("\n── Assembly ─────────────────────────────────────────────────────────")
    p23_items = list(p23_sec_amines.items())
    result_rows = []
    assemble_shared = {'p23_items': p23_items}

    # Pipeline A: P1 = secondary amine  — p23_items sent once per worker
    a_tasks = list(sec_amines_p1.items())
    print(f"  Pipeline A (P1=amine): {len(sec_amines_p1):,} × {len(p23_sec_amines):,} "
          f"= {len(sec_amines_p1)*len(p23_sec_amines):,} pairs …")
    result_rows += _parallel(_worker_assemble_amine_p1, a_tasks, workers,
                             shared=assemble_shared)
    print(f"    → {len(result_rows):,} products")

    # Pipeline B: P1 = alcohol (skippable with --skip_alcohol_p1)
    if getattr(args, 'skip_alcohol_p1', False):
        print("  Pipeline B (P1=alcohol): SKIPPED (--skip_alcohol_p1)")
    else:
        b_tasks  = list(alcohols.items())
        before_b = len(result_rows)
        print(f"  Pipeline B (P1=alcohol): {len(alcohols):,} × {len(p23_sec_amines):,} "
              f"= {len(alcohols)*len(p23_sec_amines):,} pairs …")
        result_rows += _parallel(_worker_assemble_alcohol_p1, b_tasks, workers,
                                 shared=assemble_shared)
        print(f"    → {len(result_rows)-before_b:,} products")

    result_df = pd.DataFrame(result_rows)
    print(f"\n  Total assembled: {len(result_df):,}")

    # ── 8. Save ───────────────────────────────────────────────────────────────
    out_path = os.path.join(args.enamine_dir, 'enamine_ligands.csv')
    if os.path.exists(out_path):
        out_path = os.path.join(args.enamine_dir,
                                f'enamine_ligands_{now.strftime("%Y%m%d_%H%M%S")}.csv')

    result_df[['synonyms', 'smiles']].to_csv(out_path, index=False)
    print(f"\nSaved (minimal) → {out_path}")

    full_path = out_path.replace('.csv', '_full.csv')
    result_df.to_csv(full_path, index=False)
    print(f"Saved (full)    → {full_path}")

    print("\n── Assembled molecule property summary ──────────────────────────────")
    log_pool(dict(zip(result_df['synonyms'], result_df['smiles'])),
             'assembled_molecules', out_dir)


import multiprocessing
multiprocessing.set_start_method("fork", force=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Build combinatorial triazine library from ENAMINE SDF building blocks.'
    )
    parser.add_argument('--enamine_dir', required=True,
                        help='Directory containing the SDF files.')
    # Required
    parser.add_argument('--alcohols',  required=True,
                        help='Filename prefix for alcohol/OH SDF(s).  e.g. "alcs"')
    parser.add_argument('--aldehydes', required=True,
                        help='Filename prefix for aldehyde SDF(s).  e.g. "aldehyde"')
    parser.add_argument('--amines',    required=True,
                        help='Filename prefix for amine SDF(s).  e.g. "amines"')
    # Optional
    parser.add_argument('--ester_carboxyl', default=None,
                        help='(optional) Prefix for combined ester+carboxyl SDF(s).')
    parser.add_argument('--halides', default=None,
                        help='(optional) Prefix for halide SDF(s).  Skips SN2 if absent.')
    parser.add_argument('--out_dir', default='bb_props',
                        help='Directory for per-pool property CSVs (default: bb_props/).')
    parser.add_argument('--workers', type=int, default=max(1, cpu_count() - 1),
                        help='Number of parallel worker processes '
                             f'(default: cpu_count-1 = {max(1, cpu_count()-1)}).')
    parser.add_argument('--prefilter_p23_sim', type=float, default=None,
                        help='(optional) Pre-assembly P2/P3 similarity filter: keep only '
                             'P2/P3 intermediates with Tanimoto >= this value to known TgCPL '
                             'lead P2/P3 fragments (Cpd18/12/13).  Recommended: 0.15-0.25. '
                             'Reduces the 10M+ P2/P3 pool before assembly, giving a large '
                             'speedup at the cost of some coverage.')
    parser.add_argument('--skip_alcohol_p1', action='store_true',
                        help='(optional) Skip Pipeline B (P1=alcohol × P2/P3). '
                             'Pipeline B is ~4x larger than A and can be deferred to a '
                             'second run or GreatLakes SLURM job.')
    parser.add_argument('--skip_synthesis', action='store_true',
                        help='Skip SDF loading and pool synthesis; load pre-built '
                             'pool CSVs from --out_dir and jump straight to assembly. '
                             'Requires a completed prior run up through step 6.')
    args = parser.parse_args()
    main(args)
