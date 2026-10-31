"""Legacy sparse spinor gauge and explicit export-contract regression tests."""
import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('legacy_soc', ROOT/'geometry/legacy_soc.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


class SpinorGauge(unittest.TestCase):
    def test_c3_orbit_and_fixed_atom(self):
        b = legacy.build_spinor_basis(['X']*4, {'X':'X7.0-s2p2d1f1'}, [1,2,0,3])
        u, r = b['U'], b['R']
        phases = np.repeat(np.exp(-1j*np.pi/3+2j*np.pi*np.arange(3)/3), b['dimensions'])
        self.assertLess(legacy.maxerr(r@u-u@sparse.diags(phases)), 1e-12)
        self.assertLess(legacy.maxerr(u.conj().T@u-sparse.eye(u.shape[0])), 1e-12)
        self.assertEqual(len(b['dimensions']), 3)

    def test_d3_partner_and_c2_eigen_gauge(self):
        axis = np.array([.5,np.sqrt(3)/2,0.])
        b = legacy.build_spinor_basis(['X']*4, {'X':'X7.0-s2p2d1f1'}, [1,2,0,3],
              c2perm=[0,2,1,3], q2=2*np.outer(axis,axis)-np.eye(3))
        u,c2=b['U'],b['C2'];m,n,p,q=b['dimensions']
        self.assertEqual(m,n);self.assertEqual(p,q)
        self.assertLess(legacy.maxerr(c2@u[:,:m]-u[:,m:2*m]),1e-12)
        self.assertLess(legacy.maxerr(c2@u[:,2*m:2*m+p]-1j*u[:,2*m:2*m+p]),1e-12)

    def test_nonbijective_or_bad_group_action_rejected(self):
        for perm in ([0,0,2], [1,0,2]):
            with self.assertRaises(ValueError):
                legacy.build_spinor_basis(['X']*3, {'X':'X7.0-s1'}, perm)


class MetadataContract(unittest.TestCase):
    def test_portable_d3_profile_restores_required_optimized_consumers(self):
        profile=json.loads((ROOT/'geometry/legacy_helpers/d3_fast_profile.json').read_text())
        env=profile['environment']
        for key in ('OPENMX_CPU_SOC_COMPACT','OPENMX_CPU_SOC_TR','OPENMX_CPU_SOC_REAL',
                    'OPENMX_CPU_OWNER_SPARSE','OPENMX_CPU_LOCAL_HS','OPENMX_CPU_LOCAL_OLP',
                    'OPENMX_CPU_FORCE_DM_SHARED','OPENMX_CPU_DM_BLAS3',
                    'OPENMX_CPU_DM_SORT','OPENMX_CPU_DM_TRIM','OPENMX_CPU_DM_CYCLIC','OPENMX_CPU_HNL_BATCH'):
            self.assertEqual(env[key],'1',key)
        self.assertEqual(env['OPENMX_CPU_GSOC'],'0')
        self.assertEqual(profile['soc_export_policy'],'symmetrize')
        self.assertTrue(profile['requires_metadata_environment_additions'])
        for key in ('OPENMX_CPU_GSOC_BASIS','OPENMX_CPU_C3_BASIS','OPENMX_CPU_SOC_NATIVE_MAP'):
            self.assertNotIn(key,env)
        self.assertFalse(any('/results/'in value for value in env.values()))

    def case(self, directory):
        cell=np.array([[10.,0,0],[-5,np.sqrt(3)*5,0],[0,0,20.]])
        xyz=np.array([[1,0,0],[-.5,np.sqrt(3)/2,0],[-.5,-np.sqrt(3)/2,0],[0,0,0]])
        lines=['System.Name soc','scf.SpinPolarization NC','scf.SpinOrbit.Coupling on',
               'scf.restart c2n','scf.maxIter 1','scf.EigenvalueSolver Cluster',
               'scf.Kgrid 1 1 1','MD.Type nomd','HS.fileout on','Atoms.Number 4',
               '<Definition.of.Atomic.Species','X X7.0-s1p1 X_PBE19','Definition.of.Atomic.Species>',
               'Atoms.UnitVectors.Unit Ang','<Atoms.UnitVectors',
               *[' '.join(map(str,row)) for row in cell], 'Atoms.UnitVectors>',
               'Atoms.SpeciesAndCoordinates.Unit Ang','<Atoms.SpeciesAndCoordinates',
               *[f'{i+1} X '+ ' '.join(map(str,row))+' 1.0 1.0 0.0 0.0 0.0 0' for i,row in enumerate(xyz)],
               'Atoms.SpeciesAndCoordinates>']
        inp=directory/'soc.dat';inp.write_text('\n'.join(lines)+'\n')
        q2=np.diag([1.,-1.,-1.]);ops=[]
        for q in [np.eye(3),legacy.Q3,legacy.Q3@legacy.Q3,q2,legacy.Q3@q2,legacy.Q3@legacy.Q3@q2]:
            target=xyz@q.T;perm=[int(np.argmin(np.linalg.norm(xyz-v,axis=1))) for v in target]
            w=np.rint(np.linalg.inv(cell.T)@q@cell.T).astype(int)
            ops.append({'W':w.tolist(),'t':[0.,0.,0.],'R_cart_orthogonal':q.tolist(),
                        'atom_permutation':perm,'image_shift':[[0,0,0]]*4})
        meta=directory/'symmetry.json';meta.write_text(json.dumps({'selected':{
            'valid':True,'symprec_Ang':1e-6,'map_tol_Ang':1.25e-6,'operations':ops}}))
        return inp,meta

    def test_raw_default_does_not_silently_project_exports(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);inp,meta=self.case(d)
            result=legacy.prepare(inp,meta,d/'out','D3')
            self.assertEqual(result['export_HS'],'raw')
            env=result['env_additions']
            self.assertEqual(env['OPENMX_CPU_SOC_TR'],'0')
            self.assertNotIn('OPENMX_CPU_SOC_NATIVE_GROUP',env)
            header=struct.unpack('<9i',Path(result['basis']).read_bytes()[:36])
            self.assertEqual(header[:5],(0x43334253,2,4,16,4))
            self.assertEqual(sum(header[5:]),32)

    def test_explicit_projected_d3_and_c3_map_headers(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);inp,meta=self.case(d)
            for group,ng,version in [('D3',6,2),('C3',3,1)]:
                out=d/group
                result=legacy.prepare(inp,meta,out,group,export_policy='symmetrize')
                self.assertEqual(result['export_HS'],'symmetrized')
                self.assertEqual(result['env_additions']['OPENMX_CPU_SOC_NATIVE_GROUP'],group)
                self.assertEqual(result['env_additions']['OPENMX_CPU_SOC_REAL'],'1'if group=='D3'else'0')
                self.assertEqual(result['env_additions']['OPENMX_CPU_SOC_TR'],'1'if group=='D3'else'0')
                self.assertEqual(struct.unpack('<8i',Path(result['map']).read_bytes()[:32])[:6],
                                 (0x52534433,2,4,16,ng,1))
                basis=Path(result['basis']).read_bytes()
                self.assertEqual(struct.unpack('<2i',basis[:8]),(0x43334253,version))

    def test_multistep_soc_and_existing_output_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);inp,meta=self.case(d)
            inp.write_text(inp.read_text().replace('scf.maxIter 1','scf.maxIter 2'))
            with self.assertRaisesRegex(ValueError,'one-step'):
                legacy.prepare(inp,meta,d/'out','D3')
            inp.write_text(inp.read_text().replace('scf.maxIter 2','scf.maxIter 1'))
            (d/'out').mkdir()
            with self.assertRaisesRegex(ValueError,'new'):
                legacy.prepare(inp,meta,d/'out','D3')

if __name__=='__main__':unittest.main()
