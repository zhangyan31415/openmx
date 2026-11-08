"""Regression for silent errors in the installed MUMPS tree-parallel solve.

The full material reproducer is recorded separately. This small integration
test verifies the policy at the actual EPS solve boundary and a known complex
nonorthogonal eigensystem, while retaining the parallel factorization policy.
"""
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace
import unittest


class MumpsBacksolveTest(unittest.TestCase):
    def test_safe_backsolve_and_known_nonorthogonal_spectrum(self):
        try:
            import numpy as np
            import scipy.sparse as sp
            from petsc4py import PETSc
            from slepc4py import SLEPc
        except ImportError as exc:
            self.skipTest(str(exc))
        if not PETSc.Sys.hasExternalPackage('mumps'):
            self.skipTest('MUMPS not installed')
        path=Path(os.environ.get('OPENMX_BAND_SOLVER_TEST_SOURCE',
            str(Path(__file__).resolve().parents[1]/'tools/band_solver.py')))
        spec=importlib.util.spec_from_file_location('tested_band_worker',path)
        worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
        test=self
        seen=[]

        class EPSProxy:
            def __init__(self,eps): self.eps=eps
            def __getattr__(self,key): return getattr(self.eps,key)
            def solve(self):
                self.eps.setUp()
                ksp=self.eps.getST().getKSP();ksp.setUp()
                factor=ksp.getPC().getFactorMatrix()
                seen.append(factor.getMumpsIcntl(48))
                test.assertEqual(seen[-1],0,
                    'MUMPS tree-parallel backsolve must be disabled before EPS.solve')
                return self.eps.solve()

        class EPSFactory:
            ProblemType=SLEPc.EPS.ProblemType
            Which=SLEPc.EPS.Which
            def create(self,comm): return EPSProxy(SLEPc.EPS().create(comm))

        worker.libraries=lambda:(PETSc,SimpleNamespace(EPS=EPSFactory,Sys=SLEPc.Sys),'mumps')
        n=80;rng=np.random.default_rng(3321)
        L=np.eye(n)+.03*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        U,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        exact=np.linspace(-4,3,n);S=L@L.conj().T
        H=L@((U*exact)@U.conj().T)@L.conj().T
        values,info=worker.solve(sp.csr_matrix(H),sp.csr_matrix(S),-.3,20,audit=True)
        np.testing.assert_allclose(values,exact[worker.nearest_indices(exact,-.3,20)],atol=1e-10,rtol=0)
        self.assertLess(info['residual'],1e-8)
        self.assertLess(info['Sorth'],1e-8)
        self.assertEqual(seen,[0])


if __name__=='__main__': unittest.main()
