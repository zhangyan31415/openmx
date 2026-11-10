"""Ordinary symmetry input selects qualified SOC memory paths automatically."""
import unittest
from unittest.mock import patch
from test_bootstrap import helper
class SOCMemoryDefault(unittest.TestCase):
 def setup_case(self,group='D3',route='legacy'):
  m=helper();api=m.module('soc_memory_default_api',m.CPU/'tools/input_symmetry.py')
  q=api.request({'scf.symmetry':'on'});prepared={'selected_group':group,'SOC_route':route,'NSOC_environment_additions':{},'SOC_environment_additions':{}}
  values={'scf.eigenvaluesolver':'Cluster','scf.eigen.lib':'elpa2','scf.mixing.type':'rmm-diis','scf.mixing.history':'30','scf.mixing.startpulay':'15','scf.spinpolarization':'nc','scf.spinorbit.coupling':'on','scf.restart':'c2n','scf.maxiter':'1'}
  return m,api,q,prepared,values
 def run_policy(self,group='D3',route='legacy',extra=None,threads=1,environ=None):
  m,api,q,p,v=self.setup_case(group,route);v.update(extra or {})
  with patch.dict(m.os.environ,environ or {},clear=True):env=m.profiles(p,q,api,threads=threads,values=v)['soc']
  return env,p
 def test_validated_d3_and_c3_are_on_without_environment_variables(self):
  for group in ['D3','C3']:
   env,p=self.run_policy(group)
   for key,value in [('SOC_MEMORY','1'),('VNA_FORCE_STREAM','2'),('VNA_SHARED','1'),('VNA_RAW0','1')]:self.assertEqual(env.get('OPENMX_CPU_'+key),value)
   self.assertTrue(p['SOC_memory_policy']['enabled']);self.assertEqual(p['SOC_memory_policy']['selection'],'automatic')
 def test_unqualified_paths_do_not_silently_activate(self):
  cases=[{'route':'generic'},{'group':'D2','route':'generic'},{'threads':2},
   {'extra':{'scf.mixing.type':'simple'}},{'extra':{'scf.spinpolarization':'off','scf.spinorbit.coupling':'off','scf.restart':'off','scf.maxiter':'50'}},
   {'extra':{'scf.maxiter':'2'}},{'extra':{'band.dispersion':'on'}},{'extra':{'scf.stress.tensor':'on'}},
   {'extra':{'wannier.func.calc':'on'}},{'extra':{'population.analysis.ao.wanniers':'on'}},
   {'extra':{'calc.abc.flag':'on'}},{'extra':{'cwf.fileout.flag':'on'}},{'extra':{'nho.fileout':'on'}},
   {'extra':{'scf.mixing.startpulay':'2'}}]
  for case in cases:
   with self.subTest(case=case):
    env,p=self.run_policy(**case);self.assertEqual(env.get('OPENMX_CPU_SOC_MEMORY'),'0');self.assertEqual(env.get('OPENMX_CPU_VNA_FORCE_STREAM'),'0');self.assertFalse(p['SOC_memory_policy']['enabled']);self.assertTrue(p['SOC_memory_policy']['reason'])
 def test_explicit_disable_is_retained_for_diagnostics(self):
  env,p=self.run_policy(environ={'OPENMX_CPU_SOC_MEMORY':'0'})
  self.assertEqual(env['OPENMX_CPU_SOC_MEMORY'],'0');self.assertEqual(env['OPENMX_CPU_VNA_FORCE_STREAM'],'0');self.assertEqual(p['SOC_memory_policy']['selection'],'explicit')
 def test_automatic_and_old_explicit_on_select_identical_numerical_flags(self):
  automatic,_=self.run_policy();explicit,_=self.run_policy(environ={'OPENMX_CPU_SOC_MEMORY':'1','OPENMX_CPU_VNA_FORCE_STREAM':'2'})
  self.assertEqual(automatic,explicit)
 def test_explicit_stream_off_is_respected(self):
  env,p=self.run_policy(environ={'OPENMX_CPU_VNA_FORCE_STREAM':'0'})
  self.assertEqual(env['OPENMX_CPU_SOC_MEMORY'],'1');self.assertEqual(env['OPENMX_CPU_VNA_FORCE_STREAM'],'0')
 def test_invalid_override_is_not_ignored(self):
  for setting in [{'OPENMX_CPU_SOC_MEMORY':'2'},{'OPENMX_CPU_VNA_FORCE_STREAM':'3'}]:
   with self.assertRaises(ValueError):self.run_policy(environ=setting)
 def test_bootstrap_reports_the_automatic_policy(self):
  import tempfile,json
  from pathlib import Path
  from test_usability import Usability
  m=helper()
  with tempfile.TemporaryDirectory()as tmp,patch.dict(m.os.environ,{},clear=True):
   root=Path(tmp);inp,out=Usability().fixture(root)
   result=m.bootstrap(inp,out,'soc',out/'environment.txt',runtime_cwd=root)
   self.assertTrue(result['SOC_memory_policy']['enabled'])
   self.assertEqual(result['SOC_memory_policy']['selection'],'automatic')
   report=(out/'symmetry-report.txt').read_text()
   for line in ['SOC memory policy','automatic on','VNA streaming mode','Memory policy reason']:self.assertIn(line,report)
   saved=json.loads((out/'BOOTSTRAP.json').read_text());self.assertEqual(saved['SOC_memory_policy'],result['SOC_memory_policy'])
if __name__=='__main__':unittest.main()
