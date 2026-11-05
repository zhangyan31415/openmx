"""Small physical representation checks; production construction stays orbit-local."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from scipy.linalg import eigh
P=Path(__file__).resolve().parents[1]/'geometry/prepare_soc_basis.py'
spec=importlib.util.spec_from_file_location('soc_basis',P)
m=importlib.util.module_from_spec(spec)
if P.exists(): spec.loader.exec_module(m)

class SpinorBasisTests(unittest.TestCase):
    def test_feature_available(self):
        self.assertTrue(hasattr(m,'local_sector_basis'),'general SOC basis API absent')
    def test_rotated_axis_d2(self):
        self.assertTrue(hasattr(m,'local_sector_basis'),'general SOC basis API absent')
        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
        perpendicular=np.cross(axis,[0.,0.,1.]);perpendicular/=np.linalg.norm(perpendicular)
        # Four atom regular D2 orbit, s orbitals; axes can be arbitrary.
        pz=[1,0,3,2];px=[2,3,0,1]
        dz=np.eye(4)[:,pz];dx=np.eye(4)[:,px]
        rz=np.kron(m.spin_c2(axis),dz);rx=np.kron(m.spin_c2(perpendicular),dx)
        t=m.time_reversal(4)
        q,info=m.local_sector_basis(rz,rx)
        self.assertLess(np.max(abs(q.conj().T@q-np.eye(4))),2e-12)
        self.assertLess(np.max(abs(rz@q-1j*q)),2e-12)
        self.assertLess(np.max(abs(rx@t@q.conj()-q)),2e-12)
        self.assertLess(np.max(abs(t@t.conj()+np.eye(8))),1e-15)
        allq=np.column_stack((q,t@q.conj()))
        self.assertLess(np.max(abs(allq.conj().T@allq-np.eye(8))),2e-12)
        self.assertTrue(info['reduced_real'])
    def test_nonorthogonal_generalized_problem(self):
        self.assertTrue(hasattr(m,'local_sector_basis'),'general SOC basis API absent')
        rz=np.kron(m.spin_c2([0,0,1]),np.diag([1,-1,1]))
        rx=np.kron(m.spin_c2([1,0,0]),np.diag([-1,1,1]))
        q,_=m.local_sector_basis(rz,rx);t=m.time_reversal(3)
        u=np.column_stack((q,t@q.conj()))
        a=np.array([[.2,.3,-.1],[.3,-.5,.1],[-.1,.1,.7]])
        b=np.array([[1.1,.1,.03],[.1,1.4,-.05],[.03,-.05,.9]])
        h=u@np.kron(np.eye(2),a)@u.conj().T;s=u@np.kron(np.eye(2),b)@u.conj().T
        ew,ev=eigh(h,s);er,cr=eigh(q.conj().T@h@q,q.conj().T@s@q)
        self.assertLess(np.max(abs(ew-np.sort(np.r_[er,er]))),1e-12)
        c=np.column_stack((q@cr,t@(q@cr).conj()))
        self.assertLess(np.max(abs(c.conj().T@s@c-np.eye(6))),1e-12)
        self.assertLess(np.max(abs(h@c-(s@c)*np.r_[er,er])),1e-12)
    def test_c2_complex_and_invalid_rep_fail_closed(self):
        self.assertTrue(hasattr(m,'local_sector_basis'),'general SOC basis API absent')
        r=np.kron(m.spin_c2([2,3,4]),np.eye(2));q,info=m.local_sector_basis(r)
        self.assertFalse(info['reduced_real'])
        self.assertLess(np.max(abs(r@q-1j*q)),1e-12)
        with self.assertRaises(ValueError):m.local_sector_basis(np.eye(4))
        with self.assertRaises(ValueError):m.local_sector_basis(r,r)
    def test_spinor_lift_minus_identity(self):
        self.assertTrue(hasattr(m,'spin_c2'),'general SOC basis API absent')
        for axis in [[1,0,0],[1,1,0],[-2,3,1]]:
            d=m.spin_c2(axis)
            self.assertLess(np.max(abs(d@d+np.eye(2))),1e-12)


class MetadataIntegrationTests(unittest.TestCase):
    def fixture(self,root):
        q,_=np.linalg.qr(np.array([[1.,2.,3.],[4.,1.,-2.],[2.,-3.,1.]]))
        if np.linalg.det(q)<0:q[:,0]*=-1
        cell=np.diag([4.,5.,6.])@q.T
        text='''Species.Number 1
<Definition.of.Atomic.Species
X X7.0-s1p1 X_PBE19
Definition.of.Atomic.Species>
Atoms.Number 1
Atoms.SpeciesAndCoordinates.Unit FRAC
<Atoms.SpeciesAndCoordinates
1 X 0 0 0 0.5 0.5
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
'''+''.join(' '.join(str(v) for v in row)+'\n' for row in cell)+'''Atoms.UnitVectors>
scf.SpinPolarization off
'''
        inp=root/'input.dat';inp.write_text(text)
        st=m.detector.parse_openmx(inp)
        operations=[np.eye(3,dtype=int),np.diag([-1,-1,1]),np.diag([1,-1,-1]),np.diag([-1,1,-1])]
        checked=m.detector._build_operation_set(st,operations,np.zeros((4,3)),1e-6,1.)
        self.assertTrue(checked['valid']);self.assertTrue(checked['orbital_representation']['checked'])
        import json
        meta=root/'symmetry.json';meta.write_text(json.dumps({'selected':checked}))
        return inp,meta
    def test_rotated_geometry_binary_and_subgroup(self):
        import tempfile,struct
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inp,meta=self.fixture(root)
            info=m.build(inp,meta,root/'d2.gsoc')
            self.assertEqual(info['selected_subgroup'],'D2');self.assertEqual(info['sector_sizes'],[4,4])
            self.assertEqual(info['dtype'],[1,1]);self.assertTrue(info['full_point_group_reduction'])
            with (root/'d2.gsoc').open('rb') as f:
                header=struct.unpack('<6iq',f.read(32));self.assertEqual(header[:6],(m.MAGIC,1,1,4,2,4))
                sizes=np.frombuffer(f.read(8),'<i4');dtype=np.frombuffer(f.read(8),'<i4');partner=np.frombuffer(f.read(8),'<i4')
                f.read(72+8+4*(96+4))
                cp=np.frombuffer(f.read(9*8),'<i8');ri=np.frombuffer(f.read(header[-1]*4),'<i4');cv=np.frombuffer(f.read(header[-1]*16),'<c16')
                self.assertEqual(f.read(),b'')
            from scipy.sparse import csc_matrix
            u=csc_matrix((cv,ri,cp),shape=(8,8)).toarray()
            self.assertLess(np.max(abs(u.conj().T@u-np.eye(8))),2e-12)
            self.assertLess(np.max(abs(u[:,4:]-m.time_reversal(4)@u[:,:4].conj())),2e-12)
            c2=m.build(inp,meta,root/'c2.gsoc',group='C2')
            self.assertEqual(c2['dtype'],[2,2]);self.assertFalse(c2['full_point_group_reduction'])
    def test_malformed_operations_rejected(self):
        import tempfile,json,copy
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inp,meta=self.fixture(root);valid=json.loads(meta.read_text())
            for field,value in [('atom_permutation',[1]),('image_shift',[[1,0,0]]),('W',[[2,0,0],[0,1,0],[0,0,1]])]:
                bad=copy.deepcopy(valid);bad['selected']['operations'][1][field]=value;meta.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):m.build(inp,meta,root/'bad.gsoc')
            bad=copy.deepcopy(valid);bad['selected']['operations']=bad['selected']['operations'][:-1];meta.write_text(json.dumps(bad))
            with self.assertRaises(ValueError):m.build(inp,meta,root/'bad.gsoc')

if __name__=='__main__':unittest.main()
