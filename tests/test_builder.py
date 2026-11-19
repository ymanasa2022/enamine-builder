import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from rdkit import Chem
import build_enamine as b


class BuilderTests(unittest.TestCase):
    def products(self, route, *smiles):
        return set(b.reaction_products(route, *map(b.molecule, smiles)))

    def test_reductive_amination_retains_both_ketone_branches(self):
        self.assertEqual(self.products('reductive_am', 'CN', 'O=C(c1ccc(F)cc1)C1CCCCC1'),
                         {'CNC(c1ccc(F)cc1)C1CCCCC1'})
        self.assertEqual(self.products('reductive_am', 'CN', 'O=Cc1ccccc1'), {'CNCc1ccccc1'})

    def test_core_symmetric_sites_and_nitrile_preserved(self):
        mono = self.products('snar_n', b.CORE, 'O1CCNCC1')
        self.assertEqual(len(mono), 1)
        products = self.products('snar_n', next(iter(mono)), 'CNCc1ccc(Cl)cc1')
        self.assertEqual(len(products), 1)
        m = b.molecule(next(iter(products)))
        self.assertEqual(sum(a.GetSymbol() == 'Cl' for a in m.GetAtoms()), 1)  # pendant aryl Cl untouched
        self.assertTrue(m.HasSubstructMatch(Chem.MolFromSmarts('C#N')))

    def test_transformations(self):
        self.assertEqual(self.products('sn2', 'CN', 'CCBr'), {'CCNC'})
        self.assertEqual(self.products('sn2', 'CN', 'CCF'), set())
        self.assertEqual(self.products('ester_to_aldehyde', 'CCC(=O)OC'), {'CCC=O'})
        self.assertEqual(self.products('ester_to_alcohol', 'CCC(=O)OC'), {'CCCO'})
        self.assertEqual(self.products('acid_to_aldehyde', 'CCC(=O)O'), {'CCC=O'})
        self.assertEqual(self.products('nitrile_to_amine', 'CC#N'), {'CCN'})

    def test_amides_not_classified_as_amine(self):
        self.assertFalse(b.molecule('CC(=O)NC').HasSubstructMatch(b.PAT['secondary']))
        self.assertTrue(b.molecule('CCNC').HasSubstructMatch(b.PAT['secondary']))

    def test_salts_and_stereochemistry(self):
        m = b.normalize(b.molecule('C[NH2+]C.[Cl-]'), {})
        self.assertEqual(b.canonical(m), 'CNC')
        with tempfile.TemporaryDirectory() as t:
            pool = {'a': 'CNC', 'b': 'CN(C)', 'r': 'CN[C@H](C)c1ccccc1', 's': 'CN[C@@H](C)c1ccccc1'}
            result = b.reduce_pool(pool, 'p23_arms', {}, Path(t))
            self.assertEqual(len(result), 3)

    def test_filter_validation_and_unrounded_values(self):
        with self.assertRaises(ValueError):
            b.validate_config({'inputs': {'amines': 'a.sdf'}, 'filters': {'final': {'properties': {'logD': {'max': 3}}}}})
        mw = b.properties('CNC')['MW']
        self.assertFalse(b.passes('CNC', {'properties': {'MW': {'max': mw - 0.001}}}))
        self.assertTrue(b.passes('CNC', {'properties': {'MW': None}}))
        self.assertTrue(b.passes('CNC', {'enabled': False, 'properties': {'MW': {'max': 1}}}))

    def test_diversity_reproducible(self):
        smiles = ['CNC', 'CCNC', 'O1CCNCC1', 'c1ccccc1', 'CCCCCCCC']
        a = b.diverse_indices(smiles, 3, 42)
        self.assertEqual(a, b.diverse_indices(smiles, 3, 42))
        self.assertEqual(len(set(a)), 3)

    def test_final_similarity_limit(self):
        smiles = ['CNC', 'CCNC', 'O1CCNCC1', 'c1ccccc1', 'CCCCCCCC']
        picks = b.diverse_indices(smiles, 5, 42, 0.5)
        from rdkit import DataStructs
        generator = b.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
        fps = [generator.GetFingerprint(b.molecule(smiles[i])) for i in picks]
        for i, fp in enumerate(fps):
            for other in fps[:i]:
                self.assertLessEqual(DataStructs.TanimotoSimilarity(fp, other), 0.5 + 1e-9)

    def test_secondary_alcohol_ketone_conversion(self):
        self.assertEqual(self.products('secondary_alcohol_to_ketone', 'C[C@H](O)c1ccccc1'), {'CC(=O)c1ccccc1'})
        self.assertEqual(self.products('ketone_to_secondary_alcohol', 'CC(=O)c1ccccc1'), {'CC(O)c1ccccc1'})
        self.assertEqual(self.products('secondary_alcohol_to_ketone', 'CC(C)(O)c1ccccc1'), set())
        self.assertEqual(self.products('ketone_to_secondary_alcohol', 'CC(=O)OC'), set())
        self.assertEqual(self.products('ketone_to_secondary_alcohol', 'CC(=O)NC'), set())

    def test_conversion_network_depth_cycles_and_limits(self):
        routes = {'secondary_alcohol_to_ketone': True, 'ketone_to_secondary_alcohol': True}
        expanded = b.expand_conversions({'alcohol': 'CC(O)c1ccccc1'}, routes, {'max_steps': 5})
        self.assertEqual(set(expanded.values()), {'CC(O)c1ccccc1', 'CC(=O)c1ccccc1'})
        self.assertEqual(len(expanded), 2)
        route = {'alcohol_to_aldehyde': True}
        one = b.expand_conversions({'diol': 'OCCCO'}, route, {'max_steps': 1})
        two = b.expand_conversions({'diol': 'OCCCO'}, route, {'max_steps': 2})
        self.assertNotIn('O=CCC=O', one.values())
        self.assertIn('O=CCC=O', two.values())
        with self.assertRaisesRegex(ValueError, 'Conversion limit exceeded'):
            b.expand_conversions({'diol': 'OCCCO'}, route, {'max_steps': 2, 'max_products_per_source': 1})
        self.assertEqual(b.expand_conversions({'a': 'CCO'}, route, {'max_steps': 0}), {'a': 'CCO'})

    def test_reductive_amination_preserves_other_carbonyl_types(self):
        products = self.products('reductive_am', 'CN', 'COC(=O)c1ccc(C=O)cc1')
        self.assertEqual(products, {'CNCc1ccc(C(=O)OC)cc1'})

    def test_converted_secondary_alcohol_enters_arm_building(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Functional groups, rather than source-file labels, determine roles.
            writer = Chem.SDWriter(str(root / 'mixed.sdf'))
            for smi in ['CN', 'O1CCNCC1', 'CC(O)c1ccccc1']:
                writer.write(b.molecule(smi))
            writer.close()
            c = {'inputs': {'amines': 'mixed.sdf'}, 'output_dir': 'out', 'workers': 1,
                 'routes': {'secondary_alcohol_to_ketone': True, 'p1_alcohols': False},
                 'conversions': {'max_steps': 2}}
            b.validate_config(c)
            b.run(c, root)
            with (root / 'out' / 'props_p3_carbonyls.csv').open() as f:
                carbonyls = list(csv.DictReader(f))
            self.assertEqual([r['smiles'] for r in carbonyls], ['CC(=O)c1ccccc1'])
            with (root / 'out' / 'props_p23_arms.csv').open() as f:
                arms = list(csv.DictReader(f))
            self.assertIn('CNC(C)c1ccccc1', [r['smiles'] for r in arms])
            self.assertTrue((root / 'out' / 'conversion_provenance.csv').exists())

    def test_end_to_end_serial_parallel_empty_and_existing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = Chem.SDWriter(str(root / 'amines.sdf'))
            for i, smi in enumerate(['CNC', 'C[NH2+]C.[Cl-]', 'CCNC', 'O1CCNCC1', 'CN']):
                mol = b.molecule(smi)
                mol.SetProp('CatalogID', f'id{i}')
                writer.write(mol)
            writer.close()
            writer = Chem.SDWriter(str(root / 'aldehydes.sdf'))
            writer.write(b.molecule('O=Cc1ccccc1'))
            writer.close()
            config = {'inputs': {'amines': 'amines.sdf', 'aldehydes': 'aldehydes.sdf'},
                      'routes': {'p1_alcohols': False}, 'diversity': {'final_pick_count': 3, 'final_candidate_limit': 4}}
            outputs = []
            for workers in (1, 2):
                c = copy.deepcopy(config)
                c['workers'] = workers
                c['output_dir'] = f'out{workers}'
                path = root / 'config.json'
                path.write_text(json.dumps(c))
                result = subprocess.run([sys.executable, str(Path(b.__file__).resolve()), '--config', str(path)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                with (root / f'out{workers}' / 'enamine_ligands.csv').open() as f:
                    rows = list(csv.DictReader(f))
                structures = {r['smiles'] for r in rows}
                self.assertEqual(len(rows), len(structures))
                self.assertGreater(len(rows), 0)
                outputs.append(structures)
            self.assertEqual(*outputs)
            c['output_dir'] = 'existing'
            c['diversity']['final_pick_count'] = 2
            b.select_existing(c, root, root / 'out1' / 'enamine_ligands.csv')
            with (root / 'existing' / 'enamine_ligands_selected.csv').open() as f:
                self.assertEqual(len(list(csv.DictReader(f))), 2)
            c['output_dir'] = 'empty'
            c['filters'] = {'final': {'properties': {'MW': {'max': 1}}}}
            b.run(c, root)
            with (root / 'empty' / 'enamine_ligands.csv').open() as f:
                self.assertEqual(list(csv.DictReader(f)), [])


if __name__ == '__main__':
    unittest.main()
