"""Public-option equivalence and visible runtime path selection."""
from pathlib import Path
import tempfile
import unittest
from test_bootstrap import INPUT,helper
import test_bootstrap as fixtures
class Usability(unittest.TestCase):
 def values(self):return {'scf.symmetry':'on','scf.eigenvaluesolver':'Cluster','scf.mixing.type':'rmm-diis'}
 def test_explicit_default_covariance_does_not_change_route(self):
  m=helper();values=self.values()
  self.assertEqual(m.preferred_soc_route(values),'auto')
  for value in ['1e-5','0.000010','1.0E-05']:
   self.assertEqual(m.preferred_soc_route({**values,'scf.symmetry.covariancetolerance':value}),'auto')
 def test_custom_covariance_keeps_its_contract(self):
  m=helper()
  for value in ['1e-6','1e-4','1e-10']:
   self.assertEqual(m.preferred_soc_route({**self.values(),'scf.symmetry.covariancetolerance':value}),'generic')
  advanced={k:v for k,v in self.values().items()if k!='scf.symmetry'}
  self.assertEqual(m.preferred_soc_route({**advanced,'scf.symmetry.covariancetolerance':'1e-5'}),'generic')
  self.assertEqual(m.preferred_soc_route({**advanced,'scf.symmetry.covariancetolerance':'1e-10'}),'auto')
 def test_reason_explains_generic_scope_and_custom_tolerance(self):
  m=helper()
  self.assertIn('custom covariance',m.soc_route_reason({**self.values(),'scf.symmetry.covariancetolerance':'1e-6'},'generic'))
  self.assertIn('Band.dispersion',m.soc_route_reason({**self.values(),'band.dispersion':'on'},'generic'))
  self.assertIn('specialized',m.soc_route_reason(self.values(),'legacy'))
 def fixture(self,root,covariance=None):
  text=INPUT.replace('scf.Symmetry.Mode User\n','').replace('scf.Symmetry.PointGroup D2\n','')
  text=text.replace('0 14 0','-6 10.392304845413264 0')
  text=text.replace('Atoms.Number 1','Atoms.Number 3').replace('1 H 0 0 0 .5 .5','1 H 4.2 3 7 .5 .5\n2 H 6.9 4.55884572681199 7 .5 .5\n3 H 6.9 1.44115427318801 7 .5 .5')
  text=text.replace('scf.SpinPolarization off','scf.SpinPolarization nc\nscf.SpinOrbit.Coupling on\nscf.restart c2n').replace('scf.maxIter 40','scf.maxIter 1')
  if covariance is not None:text+='scf.Symmetry.CovarianceTolerance '+covariance+'\n'
  return fixtures.Bootstrap().fixture(root,text)
 def test_real_default_route_profiles_and_summary_are_identical(self):
  results=[];m=helper()
  for covariance in [None,'1e-5']:
   with tempfile.TemporaryDirectory()as tmp:
    root=Path(tmp);inp,out=self.fixture(root,covariance)
    result=m.bootstrap(inp,out,'soc',out/'environment.txt',runtime_cwd=root)
    results.append(result)
    self.assertEqual(result['SOC_route'],'legacy')
    report=(out/'symmetry-report.txt').read_text()
    for field in ['Solver path','Path reason','Geometry tolerance','Detected at','legacy']:
     self.assertIn(field,report)
    for key in ['OPENMX_CPU_SOC_COMPACT','OPENMX_CPU_SOC_TR','OPENMX_CPU_SOC_REAL','OPENMX_CPU_OWNER_SPARSE']:
     self.assertEqual(result['environment'][key],'1')
  for key in ['OPENMX_CPU_GSOC','OPENMX_CPU_SOC_COMPACT','OPENMX_CPU_SOC_TR','OPENMX_CPU_SOC_REAL','OPENMX_CPU_OWNER_SPARSE']:
   self.assertEqual(results[0]['environment'][key],results[1]['environment'][key])
 def test_real_custom_covariance_route_is_visible(self):
  m=helper()
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);inp,out=self.fixture(root,'1e-6');result=m.bootstrap(inp,out,'soc',out/'environment.txt',runtime_cwd=root)
   self.assertEqual(result['SOC_route'],'generic')
   self.assertIn('custom covariance',result['SOC_route_reason'])
   report=(out/'symmetry-report.txt').read_text();self.assertIn('generic',report);self.assertIn('custom covariance',report)
if __name__=='__main__':unittest.main()
