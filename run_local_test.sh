#!/usr/bin/env bash
# Quick local end-to-end test for build_enamine.py.
# Creates small synthetic SDF fixtures, runs the pipeline with 2 workers,
# and verifies all output molecules contain the triazine+nitrile scaffold.
#
# Usage:  bash run_local_test.sh
#         bash run_local_test.sh /path/to/real_enamine_dir  alcs aldehyde amines halide
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TEST_DIR="${TMPDIR:-/tmp}/enamine_test_$$"
mkdir -p "$TEST_DIR"

# ── If real SDF directory supplied, use it; otherwise generate synthetic fixtures ──
if [[ $# -ge 4 ]]; then
    ENAMINE_DIR="$1"; ALC="$2"; ALD="$3"; AMI="$4"; HAL="${5:-}"; ESTER="${6:-}"
    echo "Using real SDF dir: $ENAMINE_DIR"
else
    echo "No SDF dir supplied — generating synthetic test fixtures in $TEST_DIR"
    python3 - "$TEST_DIR" << 'PY'
import sys
from rdkit import Chem

out = sys.argv[1]

def write_sdf(path, pairs):
    w = Chem.SDWriter(path)
    for i,(name,smi) in enumerate(pairs):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol.SetProp('_Name', name)
            mol.SetProp('CatalogID', f'T{i:04d}')
            w.write(mol)
    w.close()

write_sdf(f"{out}/amines_test.sdf",[
    ("methylamine","CN"),("ethylamine","CCN"),("benzylamine","NCc1ccccc1"),
    ("furfurylamine","NCc1ccco1"),("cyclohexylamine","NC1CCCCC1"),
    ("aniline","Nc1ccccc1"),("propylamine","CCCN"),("isopropylamine","CC(N)C"),
    ("allylamine","NCC=C"),("piperonylamine","NCc1ccc2c(c1)OCO2"),
])
write_sdf(f"{out}/alcs_test.sdf",[
    ("methanol","CO"),("ethanol","CCO"),("propanol","CCCO"),
    ("benzylalc","OCc1ccccc1"),("furfurylalc","OCc1ccco1"),
    ("cyclohexanol","OC1CCCCC1"),("isopropanol","CC(O)C"),("phenol","Oc1ccccc1"),
])
write_sdf(f"{out}/aldehyde_test.sdf",[
    ("benzaldehyde","O=Cc1ccccc1"),("4-fluorobenzald","O=Cc1ccc(F)cc1"),
    ("furfural","O=Cc1ccco1"),("acetaldehyde","CC=O"),("propionaldehyde","CCC=O"),
    ("isobutyraldehyde","CC(C)C=O"),("4-methoxybenzald","O=Cc1ccc(OC)cc1"),
    ("3-pyridinecarb","O=Cc1cccnc1"),("cyclohexcarb","O=CC1CCCCC1"),("pivaldehyde","CC(C)(C)C=O"),
])
write_sdf(f"{out}/halide_test.sdf",[
    ("benzyl_bromide","BrCc1ccccc1"),("allyl_bromide","BrCC=C"),
    ("methyl_iodide","CI"),("ethyl_bromide","CCBr"),
    ("3-bromopropene","BrCCC=C"),("4-fluorobenzylBr","BrCc1ccc(F)cc1"),
    ("2-bromothiophenemethyl","BrCc1cccs1"),("propargyl_bromide","BrCC#C"),
])
write_sdf(f"{out}/ester_carboxyl_test.sdf",[
    ("methyl_propionate","CCC(=O)OC"),("methyl_benzoate","COC(=O)c1ccccc1"),
    ("propionic_acid","CCC(=O)O"),("benzoic_acid","OC(=O)c1ccccc1"),
])
print("Fixtures written.")
PY
    ENAMINE_DIR="$TEST_DIR"; ALC="alcs_test"; ALD="aldehyde_test"
    AMI="amines_test"; HAL="halide_test"; ESTER="ester_carboxyl_test"
fi

OUT_DIR="$TEST_DIR/out"

echo ""
echo "▶  Running build_enamine.py with --workers 2 ..."
python3 "$SCRIPT_DIR/build_enamine.py" \
    --enamine_dir "$ENAMINE_DIR" \
    --alcohols  "$ALC" \
    --aldehydes "$ALD" \
    --amines    "$AMI" \
    ${ESTER:+--ester_carboxyl "$ESTER"} \
    ${HAL:+--halides "$HAL"} \
    --out_dir   "$OUT_DIR" \
    --workers 2

echo ""
echo "▶  Verifying scaffold integrity ..."
python3 - "$TEST_DIR/enamine_ligands.csv" << 'PY'
import sys, csv
from rdkit import Chem
tri = Chem.MolFromSmarts('c1ncncn1')
cn  = Chem.MolFromSmarts('C#N')
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f))
bad = [r['smiles'] for r in rows if not (
    (m := Chem.MolFromSmiles(r['smiles'])) and m.HasSubstructMatch(tri) and m.HasSubstructMatch(cn)
)]
print(f"Total ligands: {len(rows)}  |  Scaffold failures: {len(bad)}")
if bad:
    print("FAIL — bad SMILES:", bad[:3])
    sys.exit(1)
print("PASS ✓")
PY

echo ""
echo "Output: $TEST_DIR/enamine_ligands.csv  (and $OUT_DIR/)"
