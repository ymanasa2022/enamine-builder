"""Audit known structures against editable stage filters without enumerating."""
import argparse
import csv
import json
from pathlib import Path
import build_enamine as b


def rejection_reasons(smiles, spec):
    if not spec.get('enabled', True):
        return []
    values = b.properties(smiles)
    reasons = []
    for key, limits in spec.get('properties', {}).items():
        if limits is None:
            continue
        for side, comparison in [('min', lambda x, y: x < y), ('max', lambda x, y: x > y)]:
            bound = limits.get(side)
            if bound is not None and comparison(values[key], bound):
                reasons.append(f'{key}={values[key]:.4f} violates {side}={bound}')
    mol = b.molecule(smiles)
    for smarts in spec.get('required_smarts', []):
        if not mol.HasSubstructMatch(b.query(smarts)):
            reasons.append('missing required SMARTS: ' + smarts)
    for smarts in spec.get('excluded_smarts', []):
        if mol.HasSubstructMatch(b.query(smarts)):
            reasons.append('matched excluded SMARTS: ' + smarts)
    return reasons


def audit(config, references):
    rows = []
    for ref in references:
        stage = ref['stage']
        if stage not in b.STAGES:
            raise ValueError('Unknown reference stage: ' + stage)
        mol = b.molecule(ref['smiles'])
        if mol is None:
            raise ValueError('Invalid reference SMILES: ' + ref['name'])
        smiles = b.canonical(b.normalize(mol, config.get('standardization', {})))
        reasons = rejection_reasons(smiles, config.get('filters', {}).get(stage, {}))
        rows.append({
            'name': ref['name'], 'stage': stage, 'smiles': smiles,
            'passes': not reasons, 'rejection_reasons': '; '.join(reasons),
            'preferred': b.passes(smiles, config.get('filters', {}).get('preferred', {})) if stage == 'final' else '',
            'logP_above_4_5': b.properties(smiles)['logP'] > 4.5,
            **b.properties(smiles),
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--references', default=str(Path(__file__).with_name('sar_references.json')))
    parser.add_argument('--output', default='sar_filter_audit.csv')
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    b.validate_config(config)
    rows = audit(config, json.loads(Path(args.references).read_text()))
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['name', 'stage', 'smiles', 'passes', 'rejection_reasons', 'preferred', 'logP_above_4_5', *b.PROPERTIES])
        writer.writeheader()
        writer.writerows(rows)
    rejected = [r for r in rows if not r['passes']]
    print(f'{len(rows)} references audited; {len(rejected)} rejected. Report: {args.output}')
    for row in rejected:
        print(row['name'] + ': ' + row['rejection_reasons'])
    return bool(rejected)


if __name__ == '__main__':
    raise SystemExit(main())
