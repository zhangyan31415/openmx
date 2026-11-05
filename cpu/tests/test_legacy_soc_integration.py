"""Common preparation selects one SOC adapter and exposes export semantics."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
common=load('common_prepare_for_legacy',ROOT/'geometry/prepare.py')
fixture=load('legacy_fixture',ROOT/'tests/test_legacy_soc.py')

class CommonSOCPreparation(unittest.TestCase):
    def test_generic_default_keeps_old_basis_and_environment(self):
        def build(inp,selected,basis,metadata,group):
            Path(basis).write_bytes(b'GSOC');Path(metadata).write_text('{}')
        with tempfile.TemporaryDirectory()as d:
            d=Path(d)
            actual=common.prepare_soc_artifacts(d/'input',d/'selected',d,'D3',
                SimpleNamespace(build=build),{'OPENMX_CPU_NEIGHBOR_PAD_ANG':'0'})
            self.assertEqual(actual['route'],'generic')
            self.assertEqual(actual['export_HS'],'raw')
            self.assertEqual(actual['files']['soc_basis'],d/'spinor.gsoc')
            self.assertEqual(actual['env_additions'],{'OPENMX_CPU_GSOC':'1',
                'OPENMX_CPU_GSOC_BASIS':str(d/'spinor.gsoc'),'OPENMX_CPU_NEIGHBOR_PAD_ANG':'0'})

    def test_explicit_legacy_d3_raw_and_symmetrized_export(self):
        with tempfile.TemporaryDirectory()as d:
            d=Path(d);inp,selected=fixture.MetadataContract().case(d)
            for policy,label in [('raw','raw'),('symmetrize','symmetrized')]:
                output=d/policy;output.mkdir()
                r=common.prepare_soc_artifacts(inp,selected,output,'D3',None,{},
                    soc_route='legacy',export_policy=policy)
                self.assertEqual(r['export_HS'],label)
                self.assertEqual(r['route'],'legacy')
                self.assertTrue(r['files']['soc_metadata'].is_file())
                self.assertTrue(r['files']['soc_native_map'].is_file())
                self.assertEqual(r['env_additions']['OPENMX_CPU_SOC_TR'],'1'if policy=='symmetrize'else'0')

    def test_legacy_c3_does_not_activate_d3_antiunitary_path(self):
        with tempfile.TemporaryDirectory()as d:
            d=Path(d);inp,selected=fixture.MetadataContract().case(d)
            r=common.prepare_soc_artifacts(inp,selected,d,'C3',None,{},
                soc_route='legacy',export_policy='symmetrize')
            self.assertEqual(r['env_additions']['OPENMX_CPU_SOC_TR'],'0')
            self.assertEqual(r['env_additions']['OPENMX_CPU_SOC_REAL'],'0')
            self.assertEqual(r['env_additions']['OPENMX_CPU_D3'],'0')

    def test_explicit_generic_symmetrize_writes_native_map(self):
        with tempfile.TemporaryDirectory()as d:
            d=Path(d);inp,selected=fixture.MetadataContract().case(d)
            spinor=common.module('generic_prepare_positive',ROOT/'geometry/prepare_soc_basis.py')
            r=common.prepare_soc_artifacts(inp,selected,d,'D3',spinor,{},export_policy='symmetrize')
            self.assertEqual(r['export_HS'],'symmetrized')
            self.assertTrue(r['files']['soc_native_map'].is_file())
            self.assertEqual(r['env_additions']['OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE'],'1')
            self.assertEqual(r['env_additions']['OPENMX_CPU_GSOC_PROJECT'],'1')
            self.assertNotIn('OPENMX_CPU_SOC_NATIVE_GROUP',r['env_additions'])

    def test_unsupported_routes_groups_or_generic_projection_rejected(self):
        for group,route,policy in [('D2','legacy','raw'),('C2','legacy','symmetrize'),
                                  ('D3','unknown','raw'),
                                  ('D3','legacy','unknown')]:
            with self.subTest(group=group,route=route,policy=policy):
                with self.assertRaises(ValueError):
                    common.prepare_soc_artifacts(None,None,None,group,None,{},
                        soc_route=route,export_policy=policy)

if __name__=='__main__':unittest.main()
