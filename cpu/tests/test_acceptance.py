"""Native layout fixtures exercise complete OpenMX format12/15 fields and gates."""
from pathlib import Path
import importlib.util,json,struct,tempfile,unittest
import numpy as np
TOOLS=Path(__file__).resolve().parents[1]/'tools'
if (TOOLS/'acceptance.py').exists():
 spec=importlib.util.spec_from_file_location('acceptance',TOOLS/'acceptance.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
else:a=None

def native_fixture(directory,spin=12,imag_delta=0.):
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);path=directory/'case.scfout'
 fields=2 if spin==12 else 8;count=3 if spin==12 else 6;n=4;order=1
 data=np.zeros((fields+3*(order+1)+count,n));S=[1.,0.2,0.2,1.]
 data[0]=[-1.,0.1,0.1,-0.5];data[fields-1]=S
 if spin==15:data[1]=data[0];data[6,1]=imag_delta
 dm=data[fields+3*(order+1):];dm[0]=[.5,0.,0.,.5]
 if spin==15:dm[1]=dm[0]
 geometry=np.zeros(32);geometry[:12]=np.array([[0,8,0,0],[0,0,8,0],[0,0,0,8]]).ravel()
 with path.open('wb')as f:
  np.array([2,spin,0,0,0,0],'<i4').tofile(f);np.array([order],'<i4').tofile(f)
  np.zeros(4,'<f8').tofile(f);np.zeros(4,'<i4').tofile(f)
  np.array([1,1],'<i4').tofile(f);np.array([1,1],'<i4').tofile(f)
  np.array([1,2,2,1],'<i4').tofile(f);np.zeros(4,'<i4').tofile(f);geometry.astype('<f8').tofile(f)
  # Neighbor order atom2 is [2,1], unlike atom1 [1,2].
  reordered=data[:,[0,1,3,2]];reordered.astype('<f8').tofile(f)
 path.with_suffix('.out').write_text(' Utot. -10.0\n<coordinates.forces\n2\n1 A 0 0 0 0 0 0\n2 A 0 0 0 0 0 0\ncoordinates.forces>\n')
 return path

class NativeFormats(unittest.TestCase):
 def test_scalar_format12_all_density_fields_and_trace(self):
  with tempfile.TemporaryDirectory()as t:
   p=native_fixture(Path(t)/'a');q=native_fixture(Path(t)/'b')
   result=a.compare_outputs(p.parent,q.parent,2)
   self.assertTrue(result['passed']);self.assertEqual(result['density']['TrDS_candidate'],2)
   self.assertEqual(result['HS']['header_spin'],12)
   export=a._hs.Export(p)
   try:self.assertEqual(export.arrays[1,1],0.2)
   finally:export.close()
 def test_spinor_format15_retains_imaginary_H_and_matching_S(self):
  with tempfile.TemporaryDirectory()as t:
   p=native_fixture(Path(t)/'a',15);q=native_fixture(Path(t)/'b',15,imag_delta=0.01)
   result=a.compare_outputs(p.parent,q.parent,2)
   self.assertFalse(result['passed']);self.assertEqual(result['HS']['max_H'],0.01)
   self.assertEqual(result['HS']['max_S'],0)
 def test_truncated_metadata_rejected(self):
  with tempfile.TemporaryDirectory()as t:
   p=Path(t)/'bad.scfout';p.write_bytes(b'\x00'*8)
   with self.assertRaises(ValueError):a._hs.Export(p)
 def test_no_result_tree_code_dependency(self):
  text=(TOOLS/'acceptance.py').read_text()
  self.assertNotIn('results/openmx_',text);self.assertEqual(a.READER.parent,TOOLS)
class Gates(unittest.TestCase):
 def test_cap_normal_finish_is_not_converged(self):
  self.assertFalse(a.classify_log('The calculation was normally finished.\n<DFT> NormRD = 1e-3 Criterion = 1e-7\n<DFT> Uele = -2 dUele = 1e-3','nsoc',0)['passed'])
 def test_projected_convergence_reports_both_residuals(self):
  text='The calculation was normally finished.\n<DFT> NormRD = 3e-6 Criterion = 1e-7\n<DFT> Uele = -2 dUele = 7e-8\nCPU_NSOC_SYM_FINAL representation=PROJECTED steps=36 converged=1 dU=7e-8 NormRD_projected=9e-12 criterion=1e-7 full_steps=0\n'
  r=a.classify_log(text,'nsoc',0);self.assertTrue(r['passed']);self.assertEqual(r['convergence_representation'],'PROJECTED');self.assertEqual(r['final_NormRD'],3e-6)
 def test_missing_memory_node_cannot_create_cluster_peak(self):
  self.assertIsNone(a.summarize_memory([dict(node='a',timestamp=0,physical_bytes=10)],['a','b'])['cluster_synchronized_peak_bytes'])
 def test_peak_sum_is_not_synchronized_peak(self):
  rows=[dict(node=h,timestamp=t,physical_bytes=v)for h,t,v in [('a',1,10),('b',1,1),('a',2,1),('b',2,10)]]
  r=a.summarize_memory(rows,['a','b']);self.assertEqual(r['cluster_synchronized_peak_bytes'],11);self.assertEqual(r['sum_node_peaks_bytes'],20)
if __name__=='__main__':unittest.main()
