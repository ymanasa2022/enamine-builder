#!/usr/bin/env python3
"""Enumerate triazine nitriles using editable JSON settings. See README.md."""
import argparse
import csv
import glob
import hashlib
import heapq
import json
import math
import sqlite3
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.MolStandardize import rdMolStandardize

CORE = 'N#CC1=NC(Cl)=NC(Cl)=N1'
# Exclude amide/sulfonamide/phosphoramide nitrogens and N-heteroatom bonds.
PRIMARY = '[N;X3;H2;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-[N,O,S])]'
SECONDARY = '[N;X3;H1;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-[N,O,S])]'
ALCOHOL = '[O;X2;H1;+0]-[C;X4]'
CARBONYL = '[C;X3;H1,H0](=[O;X1])(-[#6])'
PAT = {k: Chem.MolFromSmarts(v) for k, v in {
    'primary': PRIMARY, 'secondary': SECONDARY, 'alcohol': ALCOHOL,
    'aldehyde': '[C;X3;H1](=[O;X1])-[#6]',
    'ketone': '[#6]-[C;X3;H0](=[O;X1])-[#6]',
    'halide': '[C;X4;H2,H3]-[Cl,Br,I]',
}.items()}
# Enumerate all valid sites; canonicalization removes symmetry duplicates.
TRIAZINE = '[c:1]1([Cl])[n:4][c:5][n:6][c:7][n:8]1'
RXN = {
    'snar_n': AllChem.ReactionFromSmarts(TRIAZINE + '.' + SECONDARY[:-1] + ':3]>>[c:1]1([N:3])[n:4][c:5][n:6][c:7][n:8]1'),
    'snar_o': AllChem.ReactionFromSmarts(TRIAZINE + '.[O;X2;H1:3]-[C;X4:2]>>[c:1]1([O:3][C:2])[n:4][c:5][n:6][c:7][n:8]1'),
    'reductive_am': AllChem.ReactionFromSmarts(PRIMARY[:-1] + ':1].[C;X3;H1,H0;!$(C-[N,O,S]);$([C]-[#6]):2]=[O]>>[N:1]-[C:2]'),
    'sn2': AllChem.ReactionFromSmarts(PRIMARY[:-1] + ':1].[C;X4;H2,H3:2]-[Cl,Br,I]>>[N:1]-[C:2]'),
    'alcohol_to_aldehyde': AllChem.ReactionFromSmarts('[C;X4;H2:1]-[O;H1:2]>>[CH:1]=[O:2]'),
    'secondary_alcohol_to_ketone': AllChem.ReactionFromSmarts('[C;X4;H1;$([C](-[#6])-[#6]):1]-[O;H1:2]>>[C:1]=[O:2]'),
    'ketone_to_secondary_alcohol': AllChem.ReactionFromSmarts('[C;X3;H0;$([C](-[#6])-[#6]):1]=[O:2]>>[CH:1]-[OH:2]'),
    'aldehyde_to_alcohol': AllChem.ReactionFromSmarts('[C;X3;H1;$([C]-[#6]):1]=[O:2]>>[CH2:1]-[OH:2]'),
    'ester_to_aldehyde': AllChem.ReactionFromSmarts('[#6:3]-[C:1](=[O:2])-[O]-[#6]>>[#6:3]-[CH:1]=[O:2]'),
    'acid_to_aldehyde': AllChem.ReactionFromSmarts('[#6:3]-[C:1](=[O:2])-[O;H1]>>[#6:3]-[CH:1]=[O:2]'),
    'ester_to_alcohol': AllChem.ReactionFromSmarts('[#6:3]-[C:1](=[O:2])-[O]-[#6]>>[#6:3]-[CH2:1]-[OH:2]'),
    'acid_to_alcohol': AllChem.ReactionFromSmarts('[#6:3]-[C:1](=[O:2])-[O;H1]>>[#6:3]-[CH2:1]-[OH:2]'),
    'nitrile_to_amine': AllChem.ReactionFromSmarts('[C:1]#[N:2]>>[CH2:1]-[NH2:2]'),
    'imine_to_amine': AllChem.ReactionFromSmarts('[C:1]=[N;H1:2]>>[C:1]-[NH2:2]'),
}
PROPERTIES = ('MW', 'logP', 'TPSA', 'HBD', 'HBA', 'RotBonds', 'RingCount',
              'AromaticRings', 'FormalCharge', 'FractionCSP3')
CONVERSIONS = ('alcohol_to_aldehyde', 'aldehyde_to_alcohol',
               'secondary_alcohol_to_ketone', 'ketone_to_secondary_alcohol',
               'ester_to_aldehyde', 'acid_to_aldehyde', 'ester_to_alcohol',
               'acid_to_alcohol', 'nitrile_to_amine', 'imine_to_amine')
STAGES = ('p1', 'p2_precursor', 'p3_precursor', 'p23', 'final', 'preferred')
_SHARED = {}


def molecule(smiles):
    return Chem.MolFromSmiles(smiles)


def canonical(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)


@lru_cache(maxsize=20000)
def properties(smiles):
    m = molecule(smiles)
    return {
        'MW': Descriptors.MolWt(m), 'logP': Descriptors.MolLogP(m),
        'TPSA': rdMolDescriptors.CalcTPSA(m),
        'HBD': rdMolDescriptors.CalcNumHBD(m), 'HBA': rdMolDescriptors.CalcNumHBA(m),
        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(m),
        'RingCount': rdMolDescriptors.CalcNumRings(m),
        'AromaticRings': rdMolDescriptors.CalcNumAromaticRings(m),
        'FormalCharge': Chem.GetFormalCharge(m),
        'FractionCSP3': rdMolDescriptors.CalcFractionCSP3(m),
    }


@lru_cache(maxsize=512)
def query(smarts):
    return Chem.MolFromSmarts(smarts)


def passes(smiles, spec):
    if not spec.get('enabled', True):
        return True
    p = properties(smiles)
    for key, limits in spec.get('properties', {}).items():
        if limits is None:
            continue
        lo, hi = limits.get('min'), limits.get('max')
        if lo is not None and p[key] < lo:
            return False
        if hi is not None and p[key] > hi:
            return False
    m = molecule(smiles)
    return (all(m.HasSubstructMatch(query(s)) for s in spec.get('required_smarts', []))
            and not any(m.HasSubstructMatch(query(s)) for s in spec.get('excluded_smarts', [])))


def validate_config(c):
    def keys(obj, allowed, label):
        if not isinstance(obj, dict):
            raise ValueError(f'{label} must be an object')
        extra = set(obj) - set(allowed)
        if extra:
            raise ValueError(f'Unknown {label} settings: {sorted(extra)}')
    keys(c, ('inputs', 'output_dir', 'workers', 'routes', 'standardization', 'filters', 'diversity', 'conversions'), 'top-level')
    keys(c.get('inputs', {}), ('directory', 'amines', 'alcohols', 'aldehydes', 'ketones', 'ester_carboxyl', 'halides'), 'inputs')
    if not c.get('inputs', {}).get('amines'):
        raise ValueError('inputs.amines must contain an SDF filename or glob pattern')
    for k, v in c['inputs'].items():
        if k != 'directory' and not (v is None or isinstance(v, str) or
                                      isinstance(v, list) and all(isinstance(x, str) for x in v)):
            raise ValueError(f'inputs.{k} must be a string, array of strings, or null')
    if not isinstance(c.get('workers', 1), int) or isinstance(c.get('workers'), bool) or c.get('workers', 1) < 1:
        raise ValueError('workers must be a positive integer')
    keys(c.get('routes', {}), ('p1_amines', 'p1_alcohols', 'reductive_amination', 'sn2', 'direct_secondary_amines',
         'alcohol_to_aldehyde', 'aldehyde_to_alcohol', 'ester_to_aldehyde', 'acid_to_aldehyde',
         'ester_to_alcohol', 'acid_to_alcohol', 'nitrile_to_amine', 'imine_to_amine',
         'secondary_alcohol_to_ketone', 'ketone_to_secondary_alcohol'), 'routes')
    keys(c.get('standardization', {}), ('fragment_parent', 'uncharge'), 'standardization')
    for section in ('routes', 'standardization'):
        if not all(isinstance(v, bool) for v in c.get(section, {}).values()):
            raise ValueError(f'{section} values must be true or false')
    diversity = c.get('diversity', {})
    keys(diversity, ('seed', 'pool_limits', 'final_pick_count', 'final_candidate_limit', 'final_max_similarity'), 'diversity')
    keys(diversity.get('pool_limits', {}), ('p1_amines', 'p1_alcohols', 'p2_precursors', 'p3_carbonyls', 'p3_halides', 'p23_arms'), 'diversity.pool_limits')
    for key, value in {**diversity.get('pool_limits', {}), **{k: v for k, v in diversity.items() if k not in ('pool_limits', 'final_max_similarity')}}.items():
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < (0 if key == 'seed' else 1)):
            raise ValueError(f'diversity.{key} must be a positive integer or null (seed may be zero)')
    for key in ('seed', 'final_candidate_limit'):
        if key in diversity and diversity[key] is None:
            raise ValueError(f'diversity.{key} cannot be null')
    similarity = diversity.get('final_max_similarity')
    if similarity is not None and (isinstance(similarity, bool) or not isinstance(similarity, (int, float)) or not 0 <= similarity <= 1):
        raise ValueError('final_max_similarity must be between 0 and 1, or null')
    if diversity.get('final_pick_count') and diversity.get('final_candidate_limit', 20000) < diversity['final_pick_count']:
        raise ValueError('final_candidate_limit must be at least final_pick_count')
    keys(c.get('conversions', {}), ('max_steps', 'max_products_per_source'), 'conversions')
    for key, value in c.get('conversions', {}).items():
        if isinstance(value, bool) or not isinstance(value, int) or value < (0 if key == 'max_steps' else 1):
            raise ValueError(f'conversions.{key} must be an integer within range')
    keys(c.get('filters', {}), STAGES, 'filters')
    for stage, spec in c.get('filters', {}).items():
        keys(spec, ('enabled', 'properties', 'required_smarts', 'excluded_smarts'), stage)
        if not isinstance(spec.get('enabled', True), bool):
            raise ValueError(f'{stage}.enabled must be true or false')
        keys(spec.get('properties', {}), PROPERTIES, stage + '.properties')
        for k, lim in spec.get('properties', {}).items():
            if lim is None:
                continue
            keys(lim, ('min', 'max'), stage + '.' + k)
            for value in lim.values():
                if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)):
                    raise ValueError(f'{stage}.{k}: bounds must be finite numbers or null')
            if lim.get('min') is not None and lim.get('max') is not None and lim['min'] > lim['max']:
                raise ValueError(f'{stage}.{k}: min exceeds max')
        for field in ('required_smarts', 'excluded_smarts'):
            patterns = spec.get(field, [])
            if not isinstance(patterns, list) or not all(isinstance(s, str) for s in patterns):
                raise ValueError(f'{stage}.{field} must be an array of SMARTS strings')
            for s in patterns:
                if query(s) is None:
                    raise ValueError(f'Invalid SMARTS in {stage}: {s}')


def reaction_products(name, *reactants):
    seen = set()
    for products in RXN[name].RunReactants(tuple(reactants)):
        p = products[0]
        try:
            # Explicit H counts in reaction templates are unreliable after a
            # change in bond order/degree. Recompute Hs on mapped amine/carbonyl.
            if name in ('reductive_am', 'sn2', 'snar_n', 'imine_to_amine'):
                for atom in p.GetAtoms():
                    if atom.HasProp('old_mapno') and atom.GetIntProp('old_mapno') in (1, 2, 3):
                        atom.SetNumExplicitHs(0)
                        atom.SetNoImplicit(False)
            if name == 'secondary_alcohol_to_ketone':
                for atom in p.GetAtoms():
                    if atom.HasProp('old_mapno') and atom.GetIntProp('old_mapno') in (1, 2):
                        atom.SetNumExplicitHs(0)
                        atom.SetNoImplicit(False)
            if name in ('secondary_alcohol_to_ketone', 'ketone_to_secondary_alcohol'):
                for atom in p.GetAtoms():
                    if atom.HasProp('old_mapno') and atom.GetIntProp('old_mapno') == 1:
                        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            Chem.SanitizeMol(p)
            smi = canonical(p)
        except (ValueError, RuntimeError):
            continue
        if smi not in seen:
            seen.add(smi)
            yield smi


def normalize(m, settings):
    if settings.get('fragment_parent', True):
        m = rdMolStandardize.FragmentParent(m)
    if settings.get('uncharge', True):
        m = rdMolStandardize.Uncharger().uncharge(m)
    return m


def load_inputs(c, base):
    root = Path(c['inputs'].get('directory', '.'))
    if not root.is_absolute():
        root = base / root
    pools = {}
    stats = {}
    for label in ('amines', 'alcohols', 'aldehydes', 'ketones', 'ester_carboxyl', 'halides'):
        patterns = c['inputs'].get(label) or []
        if isinstance(patterns, str):
            patterns = [patterns]
        paths = set()
        for pattern in patterns:
            matches = glob.glob(str(root / pattern))
            if not matches:
                raise FileNotFoundError(f'No files matched inputs.{label}: {root / pattern}')
            paths.update(matches)
        entries = {}
        read = invalid = 0
        for path in sorted(paths):
            for i, m in enumerate(Chem.SDMolSupplier(path, sanitize=True, removeHs=True)):
                read += 1
                if m is None:
                    invalid += 1
                    continue
                identifier = next((m.GetProp(k).strip() for k in ('CatalogID', 'ID', 'Code', '_Name')
                                   if m.HasProp(k) and m.GetProp(k).strip()), f'mol_{i}')
                try:
                    smi = canonical(normalize(m, c.get('standardization', {})))
                except (ValueError, RuntimeError):
                    invalid += 1
                    continue
                key = f'{Path(path).name}:{identifier}:{i}'
                entries[key] = smi
        pools[label] = entries
        stats[label] = {'records': read, 'invalid': invalid, 'valid': len(entries)}
        print(f'{label}: {len(entries):,} valid entries ({invalid:,} invalid)', flush=True)
    return pools, stats


def filtered(pool, c, stage):
    return {k: s for k, s in pool.items() if passes(s, c.get('filters', {}).get(stage, {}))}


def expand_conversions(source, routes, settings):
    """Bounded breadth-first conversion expansion before role classification.

    A canonical structure is visited once per source, so reversible routes do
    not loop. Keep one shortest conversion path for each source/structure.
    """
    enabled = [r for r in CONVERSIONS if routes.get(r, False)]
    max_steps = settings.get('max_steps', 2)
    max_products = settings.get('max_products_per_source', 100)
    result = {}
    for source_key, smi in source.items():
        smi = canonical(molecule(smi))
        visited = {smi}
        result[source_key] = smi
        frontier = [(source_key, smi)]
        for _ in range(max_steps):
            following = []
            for key, current in frontier:
                for route in enabled:
                    for i, product in enumerate(reaction_products(route, molecule(current))):
                        if product in visited:
                            continue
                        if len(visited) - 1 >= max_products:
                            raise ValueError(f'Conversion limit exceeded for {source_key}; '
                                             'reduce enabled routes/max_steps or increase max_products_per_source')
                        visited.add(product)
                        name = f'{key}|{route}:{i}'
                        result[name] = product
                        following.append((name, product))
            frontier = following
            if not frontier:
                break
    return result


def init_worker(shared):
    global _SHARED
    _SHARED = shared


def arm_worker(item):
    key, smi = item
    out = []
    for route, pool in _SHARED['partners']:
        for other_key, other_smi in pool:
            for i, p in enumerate(reaction_products(route, molecule(smi), molecule(other_smi))):
                if molecule(p).HasSubstructMatch(PAT['secondary']) and passes(p, _SHARED['p23_filter']):
                    out.append((f'{route}|{key}|{other_key}|site{i}', p))
    return out


def assembly_worker(task):
    p1_name, p1_smi, kind = task
    route = 'snar_n' if kind == 'amine' else 'snar_o'
    out = []
    for mono in reaction_products(route, molecule(CORE), molecule(p1_smi)):
        for p23_name, p23_smi in _SHARED['p23_items']:
            for i, full in enumerate(reaction_products('snar_n', molecule(mono), molecule(p23_smi))):
                # Ensure two substitutions on the nitrile-triazine core.
                if molecule(full).HasSubstructMatch(query('n1c(Cl)nc(C#N)nc1')):
                    continue
                if not passes(full, _SHARED['final_filter']):
                    continue
                row = {'synonyms': f'{p1_name}::{p23_name}::site{i}',
                       'p1_name': p1_name, 'p1_smiles': p1_smi, 'p1_type': kind,
                       'p23_name': p23_name, 'p23_smiles': p23_smi, 'smiles': full,
                       'preferred': passes(full, _SHARED['preferred_filter'])}
                row.update(properties(full))
                out.append(row)
    return out


def batches(worker, tasks, shared, workers):
    if workers == 1:
        init_worker(shared)
        for task in tasks:
            yield worker(task)
    else:
        with Pool(workers, initializer=init_worker, initargs=(shared,)) as pool:
            yield from pool.imap(worker, tasks, chunksize=1)


def write_pool(pool, label, out):
    with (out / f'props_{label}.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['name', 'smiles', *PROPERTIES])
        w.writeheader()
        for name, smi in pool.items():
            w.writerow({'name': name, 'smiles': smi, **properties(smi)})



def diverse_indices(smiles, count, seed, max_similarity=None):
    """Greedy MaxMin using Morgan radius-2 fingerprints; retain supplied stereo."""
    if not smiles:
        return []
    if count is None:
        count = len(smiles)
    count = min(count, len(smiles))
    if count == len(smiles) and max_similarity is None:
        return list(range(len(smiles)))
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fps = [generator.GetFingerprint(molecule(s)) for s in smiles]
    if max_similarity is not None:
        picks, _ = MaxMinPicker().LazyBitVectorPickWithThreshold(fps, len(fps), count, 1 - max_similarity, seed=seed)
        return list(picks)
    return list(MaxMinPicker().LazyBitVectorPick(fps, len(fps), count, seed=seed))


def reduce_pool(pool, label, c, out):
    # Normalize before this point, then deduplicate exact isomeric SMILES.
    unique = {}
    with (out / f'aliases_{label}.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['representative_name', 'source_name', 'smiles'])
        for name, smi in pool.items():
            smi = canonical(molecule(smi))
            representative = unique.setdefault(smi, name)
            w.writerow([representative, name, smi])
    entries = [(name, smi) for smi, name in sorted(unique.items())]
    div = c.get('diversity', {})
    count = div.get('pool_limits', {}).get(label)
    indices = diverse_indices([s for _, s in entries], count, div.get('seed', 42))
    kept = dict(entries[i] for i in indices)
    print(f'{label}: {len(pool):,} routes -> {len(unique):,} unique -> {len(kept):,} selected', flush=True)
    return kept


def select_final(db, c, out):
    div = c.get('diversity', {})
    count = div.get('final_pick_count')
    if count is None:
        return 0
    # Uniform deterministic sample of UNIQUE products using a stable hash.
    # Exact MaxMin over millions of fingerprints would be too expensive.
    limit = div.get('final_candidate_limit', 20000)
    seed = div.get('seed', 42)
    heap = []
    for (smi,) in db.execute('SELECT smiles FROM products'):
        score = int.from_bytes(hashlib.blake2b((str(seed) + ':' + smi).encode(), digest_size=16).digest(), 'big')
        item = (-score, smi)
        if len(heap) < limit:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    candidates = [s for _, s in sorted(heap, reverse=True)]
    picks = diverse_indices(candidates, count, seed, div.get('final_max_similarity'))
    with (out / 'enamine_ligands_selected.csv').open('w', newline='') as f:
        fields = ['synonyms', 'smiles', 'preferred', *PROPERTIES]
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for index in picks:
            (record,) = db.execute('SELECT record FROM products WHERE smiles=?', (candidates[index],)).fetchone()
            w.writerow(json.loads(record))
    print(f'Final selection: {len(picks):,} from {len(candidates):,} candidates', flush=True)
    return len(picks)

def run(c, base):
    out = Path(c.get('output_dir', 'triazine_output'))
    if not out.is_absolute():
        out = base / out
    out.mkdir(parents=True, exist_ok=True)
    # Avoid silently overwriting a previous enumeration.
    if (out / 'enamine_ligands.csv').exists() or (out / 'products.sqlite').exists():
        raise FileExistsError(f'{out} contains enumeration results; choose a new output_dir')
    (out / 'config_used.json').write_text(json.dumps(c, indent=2) + '\n')
    pools, stats = load_inputs(c, base)
    routes = c.get('routes', {})
    workers = c.get('workers', 1)
    # Pool labels locate files; actual functional groups determine destinations.
    source = {k: s for pool in pools.values() for k, s in pool.items()}
    variants = expand_conversions(source, routes, c.get('conversions', {}))
    stats['conversion_source_entries'] = len(source)
    stats['conversion_variants_including_sources'] = len(variants)
    with (out / 'conversion_provenance.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['source_and_conversion_path', 'smiles'])
        w.writerows(variants.items())
    carbonyls = filtered({k: s for k, s in variants.items()
                         if molecule(s).HasSubstructMatch(PAT['aldehyde']) or molecule(s).HasSubstructMatch(PAT['ketone'])}, c, 'p3_precursor')
    alcohols = filtered({k: s for k, s in variants.items() if molecule(s).HasSubstructMatch(PAT['alcohol'])}, c, 'p1') if routes.get('p1_alcohols', True) else {}
    halides = filtered({k: s for k, s in variants.items() if molecule(s).HasSubstructMatch(PAT['halide'])}, c, 'p3_precursor')
    carbonyls = reduce_pool(carbonyls, 'p3_carbonyls', c, out)
    halides = reduce_pool(halides, 'p3_halides', c, out)
    alcohols = reduce_pool(alcohols, 'p1_alcohols', c, out)
    amines = variants
    primary = filtered({k: s for k, s in amines.items() if molecule(s).HasSubstructMatch(PAT['primary'])}, c, 'p2_precursor')
    secondary = {k: s for k, s in amines.items() if molecule(s).HasSubstructMatch(PAT['secondary'])}
    p1 = filtered(secondary, c, 'p1') if routes.get('p1_amines', True) else {}
    p23 = filtered(secondary, c, 'p23') if routes.get('direct_secondary_amines', True) else {}
    primary = reduce_pool(primary, 'p2_precursors', c, out)
    p1 = reduce_pool(p1, 'p1_amines', c, out)
    partners = []
    if routes.get('reductive_amination', True):
        partners.append(('reductive_am', list(carbonyls.items())))
    if routes.get('sn2', False):
        partners.append(('sn2', list(halides.items())))
    print(f'Building P2/P3 arms from {len(primary):,} primary amines', flush=True)
    shared = {'partners': partners, 'p23_filter': c.get('filters', {}).get('p23', {})}
    for batch in batches(arm_worker, list(primary.items()), shared, workers):
        p23.update(batch)
    p23 = reduce_pool(p23, 'p23_arms', c, out)
    for label, pool in [('p1_amines', p1), ('p1_alcohols', alcohols), ('p2_precursors', primary),
                        ('p3_carbonyls', carbonyls), ('p3_halides', halides), ('p23_arms', p23)]:
        stats[label] = len(pool)
        write_pool(pool, label, out)
        print(f'{label}: {len(pool):,}', flush=True)
    tasks = [(k, s, 'amine') for k, s in p1.items()] + [(k, s, 'alcohol') for k, s in alcohols.items()]
    print(f'Assembling {len(tasks):,} P1 entries x {len(p23):,} P2/P3 arms', flush=True)
    shared = {'p23_items': list(p23.items()), 'final_filter': c.get('filters', {}).get('final', {}),
              'preferred_filter': c.get('filters', {}).get('preferred', {})}
    # Disk-backed global deduplication, including swapped P1/P23 assignments.
    db = sqlite3.connect(out / 'products.sqlite')
    db.execute('CREATE TABLE products (smiles TEXT PRIMARY KEY, record TEXT, route_count INTEGER DEFAULT 1)')
    generated = 0
    try:
        for batch in batches(assembly_worker, tasks, shared, workers):
            for row in batch:
                generated += 1
                db.execute('INSERT INTO products(smiles, record) VALUES (?, ?) ON CONFLICT(smiles) DO UPDATE SET route_count=route_count+1',
                           (row['smiles'], json.dumps(row)))
            db.commit()
        fields = ['synonyms', 'p1_name', 'p1_smiles', 'p1_type', 'p23_name', 'p23_smiles', 'smiles', 'preferred', *PROPERTIES, 'route_count']
        with (out / 'enamine_ligands.csv').open('w', newline='') as minimal, (out / 'enamine_ligands_full.csv').open('w', newline='') as full:
            w1 = csv.DictWriter(minimal, fieldnames=['synonyms', 'smiles'])
            w2 = csv.DictWriter(full, fieldnames=fields)
            w1.writeheader()
            w2.writeheader()
            for smi, record, count in db.execute('SELECT smiles, record, route_count FROM products ORDER BY smiles'):
                row = json.loads(record)
                row['route_count'] = count
                w1.writerow({'synonyms': row['synonyms'], 'smiles': smi})
                w2.writerow(row)
        stats['selected_unique'] = select_final(db, c, out)
        stats['final_routes_passing'] = generated
        stats['final_unique'] = db.execute('SELECT COUNT(*) FROM products').fetchone()[0]
    finally:
        db.close()
    (out / 'summary.json').write_text(json.dumps(stats, indent=2) + '\n')
    print(f"Saved {stats['final_unique']:,} unique products to {out}", flush=True)



def select_existing(c, base, csv_path):
    """Stream an existing CSV through normalization, filters and exact dedup."""
    out = Path(c.get('output_dir', 'triazine_output'))
    if not out.is_absolute():
        out = base / out
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'products.sqlite').exists():
        raise FileExistsError(f'{out} contains results; choose a new output_dir')
    (out / 'config_used.json').write_text(json.dumps(c, indent=2) + '\n')
    db = sqlite3.connect(out / 'products.sqlite')
    db.execute('CREATE TABLE products (smiles TEXT PRIMARY KEY, record TEXT, route_count INTEGER DEFAULT 1)')
    read = invalid = passed = 0
    try:
        with Path(csv_path).open(newline='') as f:
            reader = csv.DictReader(f)
            if 'smiles' not in (reader.fieldnames or []):
                raise ValueError('Existing CSV must have a smiles column')
            for row in reader:
                read += 1
                m = molecule(row.get('smiles', ''))
                if m is None:
                    invalid += 1
                    continue
                smi = canonical(normalize(m, c.get('standardization', {})))
                if not passes(smi, c.get('filters', {}).get('final', {})):
                    continue
                row.update({'smiles': smi, 'synonyms': row.get('synonyms') or row.get('name') or f'row_{read}',
                            'preferred': passes(smi, c.get('filters', {}).get('preferred', {})), **properties(smi)})
                db.execute('INSERT INTO products(smiles, record) VALUES (?, ?) ON CONFLICT(smiles) DO UPDATE SET route_count=route_count+1',
                           (smi, json.dumps(row)))
                passed += 1
                if passed % 10000 == 0:
                    db.commit()
                if read % 100000 == 0:
                    print(f'Read {read:,} rows; {passed:,} passing routes', flush=True)
        db.commit()
        unique = db.execute('SELECT COUNT(*) FROM products').fetchone()[0]
        picked = select_final(db, c, out)
        summary = {'input_rows': read, 'invalid': invalid, 'passing_rows': passed, 'unique_passing': unique, 'selected': picked}
        (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(f'{read:,} input rows -> {unique:,} unique passing structures -> {picked:,} selected', flush=True)
    finally:
        db.close()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True, help='JSON input file; relative paths resolve beside this file')
    p.add_argument('--select-existing', help='Select from an existing CSV instead of enumerating SDFs')
    p.add_argument('--validate-only', action='store_true', help='Validate settings without reading SDFs')
    args = p.parse_args()
    config_path = Path(args.config).resolve()
    try:
        c = json.loads(config_path.read_text())
        validate_config(c)
        if args.validate_only:
            print('Configuration valid')
        elif args.select_existing:
            select_existing(c, config_path.parent, args.select_existing)
        else:
            run(c, config_path.parent)
    except (ValueError, OSError, KeyError) as e:
        p.exit(2, f'Error: {e}\n')


if __name__ == '__main__':
    main()
