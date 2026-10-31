"""Static checks for portable profile data; no geometry solve or native job."""
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('static_profile_workflow',ROOT/'tools/workflow.py')
workflow=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(workflow)

class PortableProfileTemplates(unittest.TestCase):
    def setUp(self):
        self.document=json.loads((ROOT/'tools/examples/profile_templates.json').read_text())

    def test_four_explicit_routes_and_no_production_default(self):
        entries=self.document['templates']
        self.assertEqual(set(entries),{'D3_legacy_fast_candidate','C3_gsoc','D2_gsoc','C2_gsoc'})
        for entry in entries.values():
            self.assertFalse(entry['automatic_default'])
            self.assertFalse(entry['universal_production_certification'])
            self.assertTrue(entry['projection_requires_explicit_configuration'])

    def test_each_profile_passes_current_static_contract(self):
        for name,entry in self.document['templates'].items():
            for stage in ('nsoc','soc'):
                with self.subTest(name=name,stage=stage):
                    workflow.validate_profile(stage,entry['profiles'][stage],entry['symmetry'])

    def test_case_paths_and_automatic_covariance_relaxations_absent(self):
        for name,entry in self.document['templates'].items():
            for profile in entry['profiles'].values():
                for key,value in profile['environment'].items():
                    self.assertIsInstance(value,str)
                    self.assertNotIn('/data/',value)
                    if key.endswith(('_BASIS','_MAP')):
                        self.assertTrue(value.startswith('${PREPARED.'),(name,key,value))
                self.assertEqual(profile['environment'].get('OPENMX_CPU_NSOC_SECTOR_PARALLEL'),'0')
                self.assertNotIn('OPENMX_CPU_NSOC_SECTOR_RANKS',profile['environment'])
            soc=entry['profiles']['soc']['environment']
            if entry['soc_route']=='generic':
                self.assertEqual(soc['OPENMX_CPU_GSOC_COV_TOL'],'3e-5')
                self.assertEqual(entry['symmetry']['covariance_tolerance'],3e-5)

    def test_history_and_known_force_failures_visible(self):
        entries=self.document['templates']
        for key,history in [('D3_legacy_fast_candidate',30),('C3_gsoc',30),('D2_gsoc',35),('C2_gsoc',35)]:
            self.assertEqual(entries[key]['input_constraints']['scf.Mixing.History'],history)
        for key in ('D3_legacy_fast_candidate','C3_gsoc'):
            self.assertEqual(entries[key]['case_evidence']['raw_physical_force_gate'],'FAIL')

    def test_legacy_and_generic_compact_consumers_are_separate(self):
        for entry in self.document['templates'].values():
            env=entry['profiles']['soc']['environment']
            if entry['soc_route']=='legacy':
                self.assertEqual(env['OPENMX_CPU_SOC_TR'],'1')
                self.assertEqual(env['OPENMX_CPU_SOC_COMPACT'],'1')
                self.assertEqual(entry['symmetry']['soc_export_policy'],'symmetrize')
            else:
                self.assertEqual(env['OPENMX_CPU_SOC_TR'],'0')
                self.assertEqual(env['OPENMX_CPU_SOC_COMPACT'],'0')
                self.assertEqual(env['OPENMX_CPU_GSOC_COMPACT'],'1')
                self.assertEqual(entry['symmetry']['soc_export_policy'],'raw')

if __name__=='__main__':unittest.main()
