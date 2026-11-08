"""Sparse native H/S assembly and target selection, independent of SCF."""
import importlib.util
import os
import json
from pathlib import Path
import struct
import tempfile
import unittest
import numpy as np
import scipy.linalg

CPU = Path(__file__).resolve().parents[1]


def backend():
    path = CPU/'tools/band_solver.py'
    if not path.is_file():
        raise AssertionError('targeted band backend is missing')
    spec = importlib.util.spec_from_file_location('targeted_band_backend', path)
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


class TargetedBands(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('OPENMX_BAND_REGRESSION_DIR'), 'external saved H/S regression fixture')
    def test_saved_large_nonorthogonal_targeted_spectrum(self):
        import scipy.sparse as sp
        root=Path(os.environ['OPENMX_BAND_REGRESSION_DIR']);meta=json.loads((root/'case.json').read_text())
        m=backend();print('REGRESSION_BACKEND_LOADED',flush=True)
        H=sp.load_npz(root/'H.npz');S=sp.load_npz(root/'S.npz')
        values,audit=m.solve(H,S,meta['target_Ha'],10,audit=True)
        ref=np.load(root/'reference.npy');expected=ref[m.nearest_indices(ref,meta['target_Ha'],10)]
        self.assertEqual(len(values),10)
        # A requested-size search can return a different valid subset at a
        # degenerate boundary. Check its roots, not unrequested completeness.
        self.assertLessEqual(float(np.min(abs(values[:,None]-ref),axis=1).max()),1e-8)
        self.assertEqual(audit['search_nev'],10)
        self.assertLessEqual(audit['residual'],1e-8)
        self.assertLessEqual(audit['Sorth'],1e-8)

    def test_sparse_overlap_input_check_has_no_size_threshold_or_dense_work(self):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        from unittest.mock import patch
        m=backend()
        # Exercise both sides of the former hard-coded 16384 boundary.
        # A sparse target solve must not materialize an N x N dense check.
        for n in (41,20001):
            S=sp.eye(n,format='csr')*1.2
            with patch.object(sp.csr_matrix,'toarray',side_effect=AssertionError('unexpected dense overlap')),\
                 patch.object(spla,'eigsh',side_effect=AssertionError('unexpected extremal overlap solve')),\
                 patch.object(scipy.linalg,'cholesky',side_effect=AssertionError('unexpected dense factorization')):
                audit=m.check_overlap(S)
            self.assertEqual(audit['overlap_check_method'],'sparse input checks')
            self.assertFalse(audit['overlap_positive_definite_certified'])
            self.assertIsNone(audit['overlap_lower_bound'])
            self.assertAlmostEqual(audit['overlap_diagonal_min'],1.2)
            self.assertIsNone(audit['overlap_min'])

    def test_sparse_overlap_accepts_general_positive_basis_and_rejects_bad_diagonal(self):
        import scipy.sparse as sp
        m=backend()
        # Positive-definite localized-orbital overlap need not be diagonally dominant.
        S=sp.csr_matrix([[1.,.8,.8],[.8,1.,.8],[.8,.8,1.]])
        self.assertGreater(np.linalg.eigvalsh(S.toarray())[0],0.)
        audit=m.check_overlap(S)
        self.assertFalse(audit['overlap_positive_definite_certified'])
        self.assertIsNone(audit['overlap_lower_bound'])
        self.assertIsNone(audit['overlap_min'])
        audit=m.check_overlap(sp.eye(32,format='csr')*5e-9)
        self.assertAlmostEqual(audit['overlap_diagonal_min'],5e-9,delta=1e-16)
        for values in ([1.,1.,1e-12],[1.,1.,-.1]):
            with self.assertRaisesRegex(ValueError,'regularization'):
                m.check_overlap(sp.diags(values,format='csr'))

    def test_sparse_overlap_rejects_nonfinite_values(self):
        import scipy.sparse as sp
        m=backend()
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            m.check_overlap(sp.csr_matrix([[1.,np.nan],[np.nan,1.]]))

    def test_targeted_generalized_spectrum_after_sparse_input_checks(self):
        import scipy.sparse as sp
        m=backend();rng=np.random.default_rng(193);n=80
        L=np.tril((rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.05)+np.eye(n)
        U,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        spectrum=np.linspace(-3.,2.,n);S=L@L.conj().T
        H=L@((U*spectrum)@U.conj().T)@L.conj().T;target=-.317
        values,audit=m.solve(sp.csr_matrix(H),sp.csr_matrix(S),target,10,audit=True)
        np.testing.assert_allclose(values,spectrum[m.nearest_indices(spectrum,target,10)],atol=1e-10,rtol=0)
        self.assertLessEqual(audit['residual'],1e-8)
        self.assertLessEqual(audit['Sorth'],1e-8)
        self.assertFalse(audit['overlap_positive_definite_certified'])
        self.assertEqual(audit['search_nev'],10)
        self.assertEqual(audit['search_ncv_policy'],'min(N,max(4*nev,20))')
        self.assertEqual(audit['search_ncv'],40)

    def test_eigenpair_audit_is_opt_in_and_does_not_change_energies(self):
        import scipy.sparse as sp
        from unittest.mock import patch
        m=backend();rng=np.random.default_rng(194);n=47
        L=np.tril((rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.02)+np.eye(n)
        U,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        spectrum=np.linspace(-2.,3.,n);S=sp.csr_matrix(L@L.conj().T)
        H=sp.csr_matrix(L@((U*spectrum)@U.conj().T)@L.conj().T)
        # PETSc solves directly from CSR buffers. Python sparse products here
        # would indicate an unwanted post-solve vector audit in default mode.
        with patch.object(sp.csr_matrix,'__matmul__',side_effect=AssertionError('unexpected default eigenpair audit')):
            values,report=m.solve(H,S,.173,10)
        self.assertFalse(report['eigenpair_audit_enabled'])
        self.assertFalse(report['eigenvectors_extracted'])
        self.assertIsNone(report['residual'])
        self.assertIsNone(report['Sorth'])
        self.assertEqual(report['eigenpair_audit_seconds'],0.)
        checked,report=m.solve(H,S,.173,10,audit=True)
        np.testing.assert_allclose(values,checked,atol=1e-10,rtol=0)
        self.assertTrue(report['eigenpair_audit_enabled'])
        self.assertTrue(report['eigenvectors_extracted'])
        self.assertLessEqual(report['residual'],1e-8)
        self.assertLessEqual(report['Sorth'],1e-8)

    def test_slepc_default_reuses_csr_without_reconstructing_hermiticity(self):
        import scipy.sparse as sp
        from unittest.mock import patch
        m=backend();PETSc,_,_=m.libraries();n=47
        diagonal=np.linspace(1.,2.,n);energies=np.linspace(-2.,3.,n)
        H=sp.diags(diagonal*energies,format='csr',dtype=PETSc.ScalarType)
        S=sp.diags(diagonal,format='csr',dtype=PETSc.ScalarType)
        for a in (H,S):
            a.indices=a.indices.astype(PETSc.IntType);a.indptr=a.indptr.astype(PETSc.IntType)
        before_H=tuple(a.copy() for a in (H.indptr,H.indices,H.data))
        before_S=tuple(a.copy() for a in (S.indptr,S.indices,S.data))
        with patch.object(m,'native_hermitian',side_effect=AssertionError('unexpected solve-time matrix reconstruction')):
            values,report=m.solve(H,S,.173,10)
        np.testing.assert_allclose(values,energies[m.nearest_indices(energies,.173,10)],atol=1e-10,rtol=0)
        for current,original in zip((H.indptr,H.indices,H.data),before_H):np.testing.assert_array_equal(current,original)
        for current,original in zip((S.indptr,S.indices,S.data),before_S):np.testing.assert_array_equal(current,original)
        self.assertIsNone(report['H_hermiticity']);self.assertIsNone(report['S_hermiticity'])
        self.assertFalse(report['matrix_hermiticity_checked'])
        self.assertFalse(report['solve_time_matrix_reconstruction'])
        self.assertTrue(all(all(flags) for flags in report['csr_python_buffers_shared'].values()))

    def test_petsc_csr_buffers_share_compatible_storage_and_check_capacity(self):
        import scipy.sparse as sp
        m=backend();A=sp.eye(47,format='csr',dtype=np.complex128)
        A.indices=A.indices.astype(np.int32);A.indptr=A.indptr.astype(np.int32)
        buffers=m.petsc_csr_buffers(A,np.int32,np.complex128)
        for x,y in zip(buffers,(A.indptr,A.indices,A.data)):
            self.assertTrue(np.shares_memory(x,y));self.assertTrue(x.flags.c_contiguous)
        converted=m.petsc_csr_buffers(A,np.int64,np.complex128)
        np.testing.assert_array_equal(converted[0],A.indptr)
        self.assertEqual(converted[0].dtype,np.int64)
        with self.assertRaisesRegex(ValueError,'capacity'):
            m.petsc_csr_buffers(sp.eye(130,format='csr'),np.int8,np.complex128)

    def test_bloch_assembly_constructs_native_lower_triangle_hermitian_operator(self):
        m=backend();phase=np.exp(2j*np.pi*.173)
        data=dict(n=2,spin=0,rows=np.array([0,0,1,1]),cols=np.array([0,1,0,1]),
                  R=np.array([[0,0,0],[1,0,0],[-1,0,0],[0,0,0]]),
                  values=np.array([[1.,9.,1.5,2.],[1.,.9,.2,1.]]))
        H,S=m.bloch_matrices(data,np.array([.173,0,0]),0)
        # Preserve stock UPLO=L semantics during construction, rather than
        # rewriting an already-created CSR matrix at each solver entry.
        np.testing.assert_allclose(H.toarray(),[[1.,1.5*phase],[1.5*phase.conjugate(),2.]])
        np.testing.assert_allclose(S.toarray(),[[1.,.2*phase],[.2*phase.conjugate(),1.]])
        self.assertEqual((H-H.getH()).nnz,0);self.assertEqual((S-S.getH()).nnz,0)

    def test_full_dense_methods_on_known_nonorthogonal_spectrum(self):
        m=backend();import scipy.sparse as sp
        rng=np.random.default_rng(167);n=37
        L=np.tril((rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.02)+np.eye(n)
        U,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        expected=np.linspace(-5,3,n);S=L@L.conj().T;H=L@((U*expected)@U.conj().T)@L.conj().T
        for driver in ('gvd','evr'):
            values,audit=m.solve_dense(sp.csr_matrix(H),sp.csr_matrix(S),-.3,0,driver)
            np.testing.assert_allclose(values,expected,atol=1e-12)
            self.assertTrue(audit['full_spectrum'])
            self.assertFalse(audit['eigenvectors_computed'])
            self.assertGreater(audit['overlap_lower_bound'],1e-9)
            values,_=m.solve_dense(sp.csr_matrix(H),sp.csr_matrix(S),-.3,5,driver)
            np.testing.assert_allclose(values,expected[m.nearest_indices(expected,-.3,5)],atol=1e-12)

    def test_full_dense_rejects_native_overlap_regularization_case(self):
        m=backend();import scipy.sparse as sp
        with self.assertRaisesRegex(ValueError,'regularization'):
            m.solve_dense(sp.eye(3,format='csr'),sp.diags([1.,1.,1e-12]),0.,0,'gvd')

    def test_reference_energy_and_nearest_selection(self):
        m = backend()
        self.assertAlmostEqual(m.reference_ha('Fermi', -.2), -.2)
        self.assertAlmostEqual(m.reference_ha('-4.0', -.2), -4/27.211386245988)
        for bad in ('nan', 'inf', 'wrong'):
            with self.assertRaises(ValueError):m.reference_ha(bad, 0.)
        np.testing.assert_array_equal(m.nearest_indices(np.array([-3., -.1, .2, 1.]), 0., 2), [1, 2])

    def test_soc_bloch_assembly_and_matching_nonidentity_overlap(self):
        m = backend()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); n = 3
            # Two directed neighbor blocks, with nonzero periodic image.
            with (root/'part-000000.bin').open('wb') as stream:
                stream.write(b'OMBAND01');stream.write(struct.pack('<5i', n, 3, 2, 8, 1))
                for row,col,R in [(0,1,(1,0,0)), (1,0,(-1,0,0))]:
                    stream.write(struct.pack('<7i',row,col,1,1,*R))
                    stream.write(np.array([2.,3.,.4,.5,.1 if row==0 else -.1,.2 if row==0 else -.2,.3,.15],'<f8').tobytes())
            data = m.load_parts(root, 1)
            H,S = m.bloch_matrices(data, np.array([.125,0,0]), 0)
            phase = np.exp(2j*np.pi*.125)
            self.assertAlmostEqual(H[0,1], (2+.1j)*phase)
            self.assertAlmostEqual(H[3,4], (3+.2j)*phase)
            self.assertAlmostEqual(H[0,4], (.4+.8j)*phase)
            self.assertAlmostEqual(H[4,0], np.conj(H[0,4]))
            self.assertAlmostEqual(S[0,1], .15*phase)
            self.assertAlmostEqual(S[3,4], .15*phase)
            self.assertEqual(S[0,4], 0.)

    def test_lower_triangle_matches_native_lapack_convention(self):
        m = backend()
        import scipy.sparse as sp
        A=sp.csr_matrix(np.array([[2+1e-12j,1+2j+1e-10],[1-2j,3]]))
        result,defect=m.native_hermitian(A)
        np.testing.assert_array_equal(result.toarray(),np.array([[2,1+2j],[1-2j,3]]))
        self.assertGreater(defect,0.)
        with self.assertRaises(ValueError):m.native_hermitian(sp.csr_matrix([[1.,10.],[1.,2.]]))

    def test_worker_environment_has_no_parent_mpi_or_loader_override(self):
        source=(CPU/'overrides/nsoc/source/cpu_band_bridge.h')
        self.assertTrue(source.is_file(),'missing isolated native bridge')
        import os,shutil,subprocess,sys
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            script=root/'probe.py';script.write_text('import os\nassert not any(k.startswith(("PMI_","PMIX_","I_MPI_","OMPI_","LD_LIBRARY_PATH","LD_PRELOAD")) for k in os.environ)\nassert os.environ.get("MKL_DYNAMIC")=="TRUE", "unsafe inherited MKL dynamic policy"\nassert os.environ["MKL_NUM_THREADS"]=="4"\n')
            (root/'test.c').write_text('#include "cpu_band_bridge.h"\nint main(int n,char**v){setenv("PMI_RANK","5",1);setenv("LD_LIBRARY_PATH","/wrong/mpi",1);setenv("MKL_DYNAMIC","FALSE",1);setenv("MKL_NUM_THREADS","4",1);char*args[]={v[1],v[2],0};int rc=cpu_band_spawn(args,0);return rc||strcmp(getenv("MKL_DYNAMIC"),"FALSE");}\n')
            subprocess.run([shutil.which('cc')or'cc','-std=gnu99','-I'+str(source.parent),str(root/'test.c'),'-o',str(root/'test')],check=True,capture_output=True,text=True)
            subprocess.run([str(root/'test'),sys.executable,str(script)],check=True)

    def test_runtime_bundle_includes_standalone_band_backend(self):
        spec=importlib.util.spec_from_file_location('target_build',CPU/'build.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        with tempfile.TemporaryDirectory() as directory:
            root=m.copy_runtime_resources(Path(directory))
            self.assertTrue((root/'cpu/tools/band_solver.py').is_file())
        self.assertTrue((CPU/'environment.yml').is_file())


if __name__=='__main__':unittest.main()
