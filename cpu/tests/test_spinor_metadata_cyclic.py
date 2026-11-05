"""Cn/Dn tests exercise complete nonorthogonal spinor spectra, not sizes only."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from scipy.linalg import eigh
P=Path(__file__).resolve().parents[1]/'geometry/prepare_soc_basis.py'
spec=importlib.util.spec_from_file_location('cyclic_basis',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class CyclicPhysicalTests(unittest.TestCase):
    def test_cyclic_api_exists(self):self.assertTrue(hasattr(m,'local_cyclic_basis'),'Cn/Dn API absent')
    def physical(self,n,dihedral):
        self.assertTrue(hasattr(m,'local_cyclic_basis'),'Cn/Dn API absent')
        theta=2*np.pi/n;spin=np.diag([np.exp(-.5j*theta),np.exp(.5j*theta)])
        r=np.kron(spin,np.eye(n)[:,(np.arange(n)+1)%n])
        b=np.kron(m.spin_c2([1,0,0]),np.eye(n)[:,(-np.arange(n))%n])
        out=m.local_cyclic_basis(r,n,b if dihedral else None)
        basis=out['basis'];t=m.time_reversal(n);u=np.column_stack(basis)
        self.assertLess(m.error(u.conj().T@u-np.eye(2*n)),2e-12)
        for k,q in enumerate(basis):
            self.assertLess(m.error(r@q-out['phase'][k]*q),2e-12)
            if out['partner'][k]>=0:
                self.assertLess(m.error(q-t@basis[out['partner'][k]].conj()),2e-12)
            if abs(out['phase'][k].imag)<1e-10:
                self.assertEqual(out['partner'][k],-1,'self-conjugate sector must remain a full independent solve')
            if dihedral:self.assertLess(m.error(b@t@q.conj()-q),2e-12)
        rng=np.random.default_rng(731+n)
        def project(x):
            x=(x+x.conj().T)*.5;x=.5*(x+t@x.conj()@t.conj().T)
            matrices=[np.linalg.matrix_power(r,k) for k in range(n)]
            if dihedral:matrices += [b@g for g in matrices]
            return sum(g@x@g.conj().T for g in matrices)/len(matrices)
        z=rng.normal(size=(2*n,2*n))+1j*rng.normal(size=(2*n,2*n))
        h=project(z);s=project(z@z.conj().T)+np.eye(2*n)
        evals,vec=eigh(h,s);all_eval=[];all_vec=[];saved={}
        for k,q in enumerate(basis):
            source=out['partner'][k]
            if source<0:
                hh=q.conj().T@h@q;ss=q.conj().T@s@q
                if out['dtype'][k]==1:
                    self.assertLess(m.error(hh.imag),2e-12);self.assertLess(m.error(ss.imag),2e-12)
                    hh=hh.real;ss=ss.real
                ew,c=eigh(hh,ss);saved[k]=(ew,c)
            else:ew,c=saved[source];c=c.conj()
            all_eval.extend(ew);all_vec.append(q@c)
        c=np.column_stack(all_vec);all_eval=np.asarray(all_eval)
        self.assertLess(m.error(np.sort(all_eval)-evals),2e-12)
        self.assertLess(m.error(c.conj().T@s@c-np.eye(2*n)),2e-12)
        self.assertLess(m.error(h@c-(s@c)*all_eval),2e-12)
    def test_c3_nonorthogonal(self):self.physical(3,False)
    def test_d3_nonorthogonal(self):self.physical(3,True)
    def test_c4_nonorthogonal(self):self.physical(4,False)
    def test_d4_nonorthogonal(self):self.physical(4,True)
    def test_c6_nonorthogonal(self):self.physical(6,False)
    def test_d6_nonorthogonal(self):self.physical(6,True)


class CyclicMetadataTests(unittest.TestCase):
    def fixture(self,root,n,dihedral):
        import json
        if n in (3,6):
            cell=np.array([[5.,0,0],[-2.5,5*np.sqrt(3)/2,0],[0,0,8.]])
            w=np.array([[0,-1,0],[1,-1,0],[0,0,1]]) if n==3 else np.array([[1,-1,0],[1,0,0],[0,0,1]])
            b=np.array([[1,-1,0],[0,-1,0],[0,0,-1]])
        else:
            cell=np.diag([5.,5.,8.]);w=np.array([[0,-1,0],[1,0,0],[0,0,1]]);b=np.diag([1,-1,-1])
        text='''Species.Number 1
<Definition.of.Atomic.Species
X X7.0-s1p1d1 X_PBE19
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
        inp=root/'input.dat';inp.write_text(text);st=m.detector.parse_openmx(inp)
        ops=[np.linalg.matrix_power(w,k) for k in range(n)]
        if dihedral:ops += [b@v for v in ops]
        checked=m.detector._build_operation_set(st,ops,np.zeros((len(ops),3)),1e-6,1.)
        self.assertTrue(checked['valid']);self.assertTrue(checked['orbital_representation']['checked'])
        meta=root/'symmetry.json';meta.write_text(json.dumps({'selected':checked}));return inp,meta
    def exercise(self,n,dihedral):
        import tempfile,struct
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inp,meta=self.fixture(root,n,dihedral);group=('D' if dihedral else 'C')+str(n)
            info=m.build(inp,meta,root/'basis.gsoc',group=group)
            self.assertEqual(info['selected_subgroup'],group);self.assertEqual(sum(info['sector_sizes']),18)
            self.assertEqual(len(info['sector_sizes']),n);self.assertEqual(info['dtype'],[1 if dihedral else 2]*n)
            self.assertTrue(info['full_point_group_reduction'])
            for phase,partner in zip(info['sector_phases'],info['partner']):
                if abs(phase[1])<1e-10:self.assertEqual(partner,-1)
            auto=m.build(inp,meta,root/'auto.gsoc')
            self.assertEqual(auto['selected_subgroup'],group)
    def test_c3_metadata(self):self.exercise(3,False)
    def test_d3_metadata(self):self.exercise(3,True)
    def test_c4_metadata(self):self.exercise(4,False)
    def test_d4_metadata(self):self.exercise(4,True)
    def test_c6_metadata(self):self.exercise(6,False)
    def test_d6_metadata(self):self.exercise(6,True)

if __name__=='__main__':unittest.main()
