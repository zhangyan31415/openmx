"""Tolerance coverage and failure diagnostics; no native SCF is launched."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tempfile
import json
import unittest
from unittest.mock import patch
CPU=Path(__file__).resolve().parents[1]
def prepare_module():
 spec=importlib.util.spec_from_file_location('scan_prepare',CPU/'geometry/prepare.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
class Scan(unittest.TestCase):
 def test_scan_reaches_exact_user_limit_above_old_cap(self):
  m=prepare_module()
  for bound in (1e-3,7e-4,.01,.02,5e-5,1e-8):
   ladder=m.tolerance_ladder(bound)
   self.assertEqual(ladder[-1],bound)
   self.assertGreaterEqual(len(ladder),2)
   self.assertTrue(all(0<a<b<=bound for a,b in zip(ladder,ladder[1:])))
  self.assertIn(3e-4,m.tolerance_ladder(1e-3))
 def test_invalid_bounds_and_custom_ladders_are_rejected(self):
  m=prepare_module()
  for bound in (0,-1,float('nan'),float('inf'),.021):
   with self.assertRaises(ValueError):m.tolerance_ladder(bound)
  for custom in ([3e-5,1e-5],[1e-5,1e-5],[float('nan')],[0]):
   with self.assertRaises(ValueError):m.tolerance_ladder(1e-3,custom)
 def test_custom_scan_includes_user_endpoint_and_upper_platform_probe(self):
  m=prepare_module();ladder=m.tolerance_ladder(.001,[1e-6,1e-5]);self.assertEqual(ladder[-1],.001);self.assertTrue(any(.0001<x<.001 for x in ladder))
 def record(self,tol,number=1,valid=True):
  return dict(valid=valid,symprec_Ang=tol,operation_count=number,max_geometry_residual_Ang=tol/4,point_group={'symbol':'32'if number>1 else'1'},operations=[{}]*number)
 def detector(self,fn):
  return SimpleNamespace(parse_openmx=lambda p:{'fake':True},detect_at_tolerance=lambda st,t,mf:fn(t))
 def test_looser_candidate_is_diagnostic_only(self):
  m=prepare_module();initial=self.record(1e-4);full={'selected':initial,'scan':[initial]};d=self.detector(lambda t:self.record(t,6 if t>=3e-4 else 1))
  def select(record,group,spinor):
   if record['operation_count']==1:raise ValueError('identity only')
   return 'D3',record['operations']
  with patch.object(m,'select_subgroup',side_effect=select):
   report=m.failure_diagnostic(Path('input.dat'),full,[1e-5,1e-4],1e-4,d,None,'auto')
  self.assertEqual(report['classification'],'larger_tolerance_candidate');self.assertEqual(report['requested_max_Ang'],1e-4);self.assertGreater(report['candidate']['symprec_Ang'],1e-4);self.assertFalse(report['runtime_launch_allowed']);self.assertFalse(report['geometry_modified']);self.assertTrue(initial==full['selected'])
 def test_identity_only_does_not_claim_proven_absence(self):
  m=prepare_module();initial=self.record(1e-4);d=self.detector(lambda t:self.record(t))
  with patch.object(m,'select_subgroup',side_effect=ValueError('identity only')):
   report=m.failure_diagnostic(Path('input.dat'),{'selected':initial,'scan':[initial]},[1e-5,1e-4],1e-4,d,None,'auto')
  self.assertEqual(report['classification'],'identity_only_in_scanned_range');self.assertFalse(report['absence_proven']);self.assertIn('scf.Symmetry off',m.diagnostic_message(report))
 def test_validation_failure_is_not_labelled_no_symmetry(self):
  m=prepare_module();bad=self.record(1e-4,6,False);bad['error']='atom map is not a bijection';d=self.detector(lambda t:self.record(t,6,False))
  report=m.failure_diagnostic(Path('input.dat'),{'selected':bad,'scan':[bad]},[1e-5,1e-4],1e-4,d,None,'auto')
  self.assertEqual(report['classification'],'operation_validation_failed');self.assertFalse(report['runtime_launch_allowed'])
 def test_unsupported_group_is_not_labelled_identity(self):
  m=prepare_module();checked=self.record(1e-4,6);d=self.detector(lambda t:checked)
  with patch.object(m,'select_subgroup',side_effect=ValueError('unsupported group')):
   report=m.failure_diagnostic(Path('input.dat'),{'selected':checked,'scan':[checked]},[1e-5,1e-4],1e-4,d,None,'auto')
  self.assertEqual(report['classification'],'unsupported_subgroup')
 def test_diagnostic_probes_stop_at_absolute_cap(self):
  m=prepare_module();rec=self.record(.02);called=[]
  d=self.detector(lambda t:called.append(t)or self.record(t))
  report=m.failure_diagnostic(Path('input.dat'),{'selected':rec,'scan':[rec]},[.01,.02],.02,d,None,'auto')
  self.assertEqual(called,[]);self.assertEqual(report['diagnostic_scan'],[]);self.assertFalse(report['runtime_launch_allowed'])
 def test_probe_errors_remain_visible_without_suggesting_a_tolerance(self):
  m=prepare_module();rec=self.record(1e-4);d=self.detector(lambda t:(_ for _ in ()).throw(ValueError('bad map')))
  report=m.failure_diagnostic(Path('input.dat'),{'selected':rec,'scan':[rec]},[1e-5,1e-4],1e-4,d,None,'auto')
  self.assertTrue(all(not r['valid'] and r['error']=='bad map' for r in report['diagnostic_scan']));self.assertIsNone(report['candidate'])
class RealGeometry(unittest.TestCase):
 def bootstrap_fixture(self, root, cell, positions, tolerance):
  from test_bootstrap import INPUT,helper,Bootstrap
  text=INPUT.replace('scf.Symmetry.Mode User\n','').replace('scf.Symmetry.PointGroup D2\n','')
  text=text.replace('12 0 0\n0 14 0\n0 0 16',cell).replace('Atoms.Number 1','Atoms.Number 3')
  text=text.replace('1 H 0 0 0 .5 .5',positions)+f'\nscf.Symmetry.Tolerance {tolerance}\n'
  return helper(),Bootstrap().fixture(root,text)
 def test_real_P1_failure_has_actionable_report_and_no_launch_marker(self):
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);m,(inp,out)=self.bootstrap_fixture(root,'11 0 0\n1.2 13 0\n.7 1.1 17','1 H .23 1.17 3.25 .5 .5\n2 H 3.81 2.47 6.12 .5 .5\n3 H 7.41 8.33 10.91 .5 .5',1e-4)
   original=inp.read_bytes()
   with self.assertRaisesRegex(ValueError,'identity_only_in_scanned_range'):
    m.bootstrap(inp,out,'nsoc',out/'environment.txt',runtime_cwd=root)
   report=json.loads((out/'metadata/SYMMETRY_DIAGNOSTIC.json').read_text())
   self.assertFalse(report['runtime_launch_allowed']);self.assertFalse((out/'environment.txt').exists());self.assertEqual(original,inp.read_bytes())
 def test_real_perturbed_D3_needs_above_old_1e4_cap(self):
  positions='1 H 4.2006 3.0001 7.0004 .5 .5\n2 H 6.8997 4.5592457268 6.9998 .5 .5\n3 H 6.9002 1.4406542732 7.0001 .5 .5'
  cell='12 0 0\n-6 10.392304845413264 0\n0 0 16'
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);m,(inp,out)=self.bootstrap_fixture(root,cell,positions,1e-4)
   with self.assertRaisesRegex(ValueError,'larger_tolerance_candidate'):
    m.bootstrap(inp,out,'nsoc',out/'environment.txt',runtime_cwd=root)
   report=json.loads((out/'metadata/SYMMETRY_DIAGNOSTIC.json').read_text())
   self.assertGreater(report['candidate']['symprec_Ang'],1e-4)
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);m,(inp,out)=self.bootstrap_fixture(root,cell,positions,3e-3);original=inp.read_bytes()
   result=m.bootstrap(inp,out,'nsoc',out/'environment.txt',runtime_cwd=root)
   self.assertEqual(result['selected_group'],'D3');self.assertEqual(original,inp.read_bytes())
   prepared=json.loads((out/'metadata/PREPARED.json').read_text());self.assertEqual(prepared['solve_tolerance_ladder_Ang'][-1],3e-3)

if __name__=='__main__':unittest.main()
