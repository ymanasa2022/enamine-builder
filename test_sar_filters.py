import copy
import json
from pathlib import Path
import unittest
import check_sar_filters as audit
import build_enamine as b


class SarTests(unittest.TestCase):
    def setUp(self):
        root = Path(__file__).parent
        self.config = json.loads((root / 'filters.json').read_text())
        self.refs = json.loads((root / 'sar_references.json').read_text())

    def test_reference_stages_pass(self):
        b.validate_config(self.config)
        for row in audit.audit(self.config, self.refs):
            self.assertTrue(row['passes'], row)

    def test_old_limits_reject_selective_references(self):
        old = copy.deepcopy(self.config)
        old['filters']['final']['properties']['logP']['max'] = 4.5
        old['filters']['p23']['properties']['RotBonds']['max'] = 7
        rows = {r['name']: r for r in audit.audit(old, self.refs)}
        self.assertIn('logP', rows['13 final']['rejection_reasons'])
        self.assertIn('RotBonds', rows['12 free arm']['rejection_reasons'])

    def test_preferred_is_not_hard_exclusion(self):
        rows = {r['name']: r for r in audit.audit(self.config, self.refs)}
        self.assertTrue(rows['18 final']['passes'])
        self.assertFalse(rows['18 final']['preferred'])
        self.assertTrue(rows['13 final']['logP_above_4_5'])

    def test_tfa_salt_normalizes_to_reference(self):
        ref = next(r for r in self.refs if r['name'] == '18 final')
        salt = dict(ref, smiles=ref['smiles'] + '.O=C(O)C(F)(F)F')
        rows = audit.audit(self.config, [ref, salt])
        self.assertEqual(rows[0]['smiles'], rows[1]['smiles'])
        self.assertEqual(rows[0]['MW'], rows[1]['MW'])


if __name__ == '__main__':
    unittest.main()
