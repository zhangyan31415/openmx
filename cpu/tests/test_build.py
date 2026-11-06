"""Meaningful standalone-source build contract checks (no native compilation)."""
import importlib.util,json,re,subprocess,sys,tempfile,unittest
from pathlib import Path
SPEC=importlib.util.spec_from_file_location('cpu_build',Path(__file__).resolve().parents[1]/'build.py')
B=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(B)
class BuildContract(unittest.TestCase):
 def test_cli_default_is_complete_unified_dry_run(self):
  with tempfile.TemporaryDirectory(prefix='openmx-build-cli-')as directory:
   output=Path(directory)/'new-output'
   run=subprocess.run([sys.executable,str(B.REPO/'cpu/build.py'),'--output',str(output)],capture_output=True,text=True,check=True)
   plan=json.loads(run.stdout)
   self.assertEqual(plan['profile'],'unified');self.assertEqual(plan['state'],'planned')
   self.assertEqual(len(plan['compile_commands']),311)
   self.assertTrue(plan['capabilities']['nsoc_symmetry']);self.assertTrue(plan['capabilities']['generic_soc'])
   self.assertEqual(Path(plan['selection']),B.REPO/'cpu/overrides/nsoc/UNIFIED_SELECTION.json')
   self.assertEqual(plan['external_openmx_objects'],[]);self.assertEqual(plan['external_elpa_objects'],[])
   self.assertFalse(output.exists())
 def test_all_sources_are_repo_local_and_text(self):
  plan=B.source_plan()
  self.assertEqual(len(plan['units']),300)
  for u in plan['units']:
   self.assertTrue((B.REPO/u['source']).is_file());self.assertNotIn('results/',u['source']);self.assertNotIn('..',Path(u['source']).parts)
 def test_stock_profile_has_no_optimized_objects(self):
  plan=B.source_plan(profile='stock')
  self.assertTrue(all('/vendor/'in u['source']for u in plan['units']))
  self.assertIn('Cluster_DFT_NonCol.o',[u['object']for u in plan['units']]);self.assertIn('DFT.o',[u['object']for u in plan['units']]);self.assertFalse(plan['capabilities']['generic_soc'])
 def test_stock_profile_includes_stock_make_targets(self):
  text=(B.REPO/'cpu/vendor/openmx4.0/source/makefile').read_text()
  objs=re.search(r'^OBJS\s*=((?:.*\\\n)*.*)',text,re.M).group(1)
  expected=set(re.findall(r'[A-Za-z0-9_]+\.o',objs))
  actual={u['object']for u in B.source_plan(profile='stock')['units']}
  self.assertFalse(expected-actual)
 def test_unified_scalar_capability_from_selected_helper(self):
  selection=json.loads((B.REPO/'cpu/overrides/nsoc/UNIFIED_SELECTION.json').read_text())
  plan=B.source_plan(selection,profile='unified')
  self.assertTrue(plan['capabilities']['nsoc_symmetry'])
  merger=next(u for u in plan['units']if u['object']=='Merge_LogFile.o')
  self.assertEqual(merger['source'],'cpu/overrides/nsoc/source/Merge_LogFile.symmetry.c')
 def test_no_binary_vendor_payload(self):
  forbidden={'.o','.mod','.a','.so'}
  self.assertFalse([p for p in (B.REPO/'cpu/vendor').rglob('*')if p.suffix in forbidden])
 def test_selection_rejects_external_source(self):
  with self.assertRaises(ValueError): B.source_plan({'replacements':{'DFT.compact.o':'/external/DFT.c'}})
 def test_selection_rejects_unknown_unit(self):
  with self.assertRaises(ValueError): B.source_plan({'replacements':{'typo.o':'cpu/overrides/soc/source/DFT.compact.c'}})
 def test_kernel_keeps_provider_dispatch(self):
  unit=next(u for u in B.source_plan()['units']if u['object']=='cpu_dm_compact_trim.gsoc.o')
  self.assertIn('-DGSOC_DENSITY_REMAP_C3_API',unit['flags']);self.assertIn('-include',unit['flags'])
if __name__=='__main__':unittest.main()
