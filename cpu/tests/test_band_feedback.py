"""Regression for roundoff overlap diagonals and masked SLEPc errors."""
import sys
from pathlib import Path
import importlib.util
import unittest
import numpy as np
import scipy.sparse as sp
CPU=Path(__file__).resolve().parents[1]
def backend():
 spec=importlib.util.spec_from_file_location('feedback_band',CPU/'tools/band_solver.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
class BandFeedback(unittest.TestCase):
 def test_roundoff_imaginary_diagonal_is_not_a_physical_failure(self):
  m=backend();S=sp.eye(4,dtype=complex,format='csr');S[2,2]+=1.43e-18j
  before=S.copy();audit=m.check_overlap(S);self.assertLess(audit['overlap_diagonal_imag_max'],1e-15);self.assertEqual((S-before).nnz,0)
 def test_real_diagonal_defect_is_still_rejected(self):
  m=backend();S=sp.eye(4,dtype=complex,format='csr');S[1,1]+=1e-8j
  with self.assertRaisesRegex(ValueError,'diagonal is not real'):m.check_overlap(S)
 def test_stock_diamond_gamma_targeted_bands_match_generalized_reference(self):
  # Native H/S from stock work/Cdia.dat. Exactly Hermitian, S positive definite;
  # the old 10*eps BV guard rejected ~1.7e-14 imaginary roundoff at Gamma.
  import scipy.linalg as la
  m=backend()
  try:m.libraries()
  except ImportError:self.skipTest('complex PETSc/SLEPc not installed')
  with np.load(CPU/'tests/data/cdia_gamma.npz',allow_pickle=False)as saved:
   H=sp.csr_matrix(saved['H']);S=sp.csr_matrix(saved['S']);target=float(saved['target'])
  before_H=H.copy();before_S=S.copy()
  ref=la.eigvalsh(H.toarray(),S.toarray());expected=ref[m.nearest_indices(ref,target,4)]
  values,report=m.solve(H,S,target,4,audit=True)
  np.testing.assert_allclose(values,expected,atol=1e-10,rtol=0)
  self.assertLess(report['residual'],1e-8);self.assertLess(report['Sorth'],1e-8)
  self.assertEqual(report['requested'],4);self.assertEqual(report['search_nev'],4)
  self.assertLess(report['basis_inner_product_tolerance'],1e-12)
  self.assertEqual((H-before_H).nnz,0);self.assertEqual((S-before_S).nnz,0)
 def test_cleanup_does_not_mask_original_error_and_cleans_every_object(self):
  m=backend();called=[]
  class Bad:
   def destroy(self):called.append('bad');raise RuntimeError('destroy failed')
  class Good:
   def destroy(self):called.append('good')
  original=ValueError('Krylov-Schur orthogonalization failed')
  with self.assertRaisesRegex(ValueError,'Krylov-Schur orthogonalization')as ctx:
   try:raise original
   finally:m.destroy_petsc_objects([Bad(),Good(),Bad(),None],sys.exc_info()[1])
  self.assertIs(ctx.exception,original);self.assertEqual(called,['bad','good','bad']);self.assertTrue(any('destroy failed'in n for n in original.__notes__))
 def test_cleanup_failure_without_solver_error_is_not_hidden(self):
  m=backend()
  class Bad:
   def destroy(self):raise RuntimeError('cleanup only failure')
  with self.assertRaisesRegex(RuntimeError,'cleanup only failure'):m.destroy_petsc_objects([Bad()],None)
if __name__=='__main__':unittest.main()
