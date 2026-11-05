"""Improper spatial operations retain polar AO and axial spin conventions."""
import importlib.util
from pathlib import Path
import unittest,numpy as np
from scipy.linalg import eigh
P=Path(__file__).resolve().parents[1]/'geometry/prepare_soc_basis.py';spec=importlib.util.spec_from_file_location('improper_basis',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ImproperTests(unittest.TestCase):
    def physical(self,kind):
        self.assertTrue(hasattr(m,'local_improper_basis'),'Ci/Cs API absent')
        if kind=='Ci':r=np.kron(np.eye(2),np.diag([1,-1,1,-1]))
        else:
            normal=np.array([1.,2.,3.]);normal/=np.linalg.norm(normal);mirror=np.eye(3)-2*np.outer(normal,normal)
            spin=m.spin_rotation(-mirror);r=np.kron(spin,np.diag([1,-1,1,-1]))
            self.assertLess(m.error(spin@spin+np.eye(2)),1e-12)
        out=m.local_improper_basis(r,kind);basis=out['basis'];u=np.column_stack(basis);t=m.time_reversal(4)
        self.assertLess(m.error(u.conj().T@u-np.eye(8)),1e-12)
        self.assertEqual(out['partner'],[-1,-1] if kind=='Ci' else [-1,0])
        self.assertEqual(out['dtype'],[2,2])
        rng=np.random.default_rng(729)
        def project(x):
            x=(x+x.conj().T)*.5;x=.5*(x+t@x.conj()@t.conj().T)
            return .5*(x+r@x@r.conj().T)
        z=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));h=project(z);s=project(z@z.conj().T)+np.eye(8)
        full=eigh(h,s,eigvals_only=True);es=[];vs=[];saved={}
        for k,q in enumerate(basis):
            if out['partner'][k]<0:
                e,c=eigh(q.conj().T@h@q,q.conj().T@s@q);saved[k]=(e,c)
            else:e,c=saved[out['partner'][k]];c=c.conj()
            es.extend(e);vs.append(q@c)
        es=np.array(es);c=np.column_stack(vs)
        self.assertLess(m.error(np.sort(es)-full),1e-12)
        self.assertLess(m.error(h@c-(s@c)*es),1e-12)
        self.assertLess(m.error(c.conj().T@s@c-np.eye(8)),1e-12)
    def test_ci_full_parity_nonorthogonal(self):self.physical('Ci')
    def test_cs_axial_spin_nonorthogonal(self):self.physical('Cs')
    def test_improper_metadata_v2(self):
        self.assertTrue(hasattr(m,'local_improper_basis'),'Ci/Cs API absent')
        import tempfile,json,struct
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inp=root/'input.dat'
            inp.write_text('''Species.Number 1
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
4 0 0
0 5 0
0 0 6
Atoms.UnitVectors>
scf.SpinPolarization off
''')
            st=m.detector.parse_openmx(inp)
            for kind,op in [('Ci',-np.eye(3,dtype=int)),('Cs',np.diag([1,1,-1]))]:
                selected=m.detector._build_operation_set(st,[np.eye(3,dtype=int),op],np.zeros((2,3)),1e-6,1.)
                self.assertTrue(selected['valid']);meta=root/f'{kind}.json';meta.write_text(json.dumps({'selected':selected}))
                out=root/f'{kind}.gsoc';info=m.build(inp,meta,out,group=kind)
                self.assertEqual(info['selected_subgroup'],kind);self.assertEqual(info['binary_version'],2)
                with out.open('rb') as f:self.assertEqual(struct.unpack('<2i',f.read(8)),(m.MAGIC,2))
                self.assertEqual(info['partner'],[-1,-1] if kind=='Ci' else [-1,0])
                self.assertEqual(sum(info['sector_sizes']),18)
                auto=m.build(inp,meta,root/f'{kind}_auto.gsoc')
                self.assertEqual(auto['selected_subgroup'],kind)
if __name__=='__main__':unittest.main()
