"""Sparse native H/S assembly and target selection, independent of SCF."""
import importlib.util
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
            script=root/'probe.py';script.write_text('import os\nassert not any(k.startswith(("PMI_","PMIX_","I_MPI_","OMPI_","LD_LIBRARY_PATH","LD_PRELOAD")) for k in os.environ)\n')
            (root/'test.c').write_text('#include "cpu_band_bridge.h"\nint main(int n,char**v){setenv("PMI_RANK","5",1);setenv("LD_LIBRARY_PATH","/wrong/mpi",1);char*args[]={v[1],v[2],0};return cpu_band_spawn(args,0);}\n')
            subprocess.run([shutil.which('cc')or'cc','-std=gnu99','-I'+str(source.parent),str(root/'test.c'),'-o',str(root/'test')],check=True,capture_output=True,text=True)
            subprocess.run([str(root/'test'),sys.executable,str(script)],check=True)

    def test_runtime_bundle_includes_standalone_band_backend(self):
        spec=importlib.util.spec_from_file_location('target_build',CPU/'build.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        with tempfile.TemporaryDirectory() as directory:
            root=m.copy_runtime_resources(Path(directory))
            self.assertTrue((root/'cpu/tools/band_solver.py').is_file())
        self.assertTrue((CPU/'environment.yml').is_file())


if __name__=='__main__':unittest.main()
