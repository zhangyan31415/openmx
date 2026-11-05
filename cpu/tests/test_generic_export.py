"""Native generic spinor export contracts and metadata regressions."""
import importlib.util,json,struct,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
class GenericExport(unittest.TestCase):
 def test_explicit_generic_symmetrize_is_supported(self):
  p=load('generic_export_prepare',ROOT/'geometry/prepare.py')
  p._soc_options('D2','generic','symmetrize')
 def test_axial_spin_lift_and_more_than_six_operations(self):
  helper=load('generic_export_map',ROOT/'geometry/generic_soc_export.py')
  rotations=[]
  for k in range(8):
   a=k*np.pi/4;rotations.append(np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]]))
  rotations.extend([-q for q in rotations.copy()])
  lifts=helper.spin_lifts(rotations)
  pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.])]
  for q,u in zip(rotations,lifts):
   for j,s in enumerate(pauli):
    expected=sum(np.linalg.det(q)*q[i,j]*pauli[i]for i in range(3))
    self.assertLess(np.max(abs(u@s@u.conj().T-expected)),1e-12)
  self.assertEqual(len(lifts),16)
 def test_group_average_matches_solved_sector_TR_real_operator(self):
  helper=load('generic_export_equivalence',ROOT/'geometry/generic_soc_export.py')
  spinor=helper.spinor;rng=np.random.default_rng(42);j=spinor.time_reversal(4)
  raw=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));raw=(raw+raw.conj().T)*.5
  raw=(raw+j@raw.conj()@j.conj().T)*.5
  for order,dihedral in [(2,False),(2,True),(3,False),(3,True)]:
   with self.subTest(order=order,dihedral=dihedral):
    angle=2*np.pi/order;q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    q2=np.diag([1.,-1.,-1.])
    def representation(q):
     orbital=spinor.detector._orbital_rotation_block('X7.0-s1p1',q)
     return np.kron(spinor.spin_rotation(q),orbital)
    rep=representation(q);perp=representation(q2)if dihedral else None
    layout=spinor.local_cyclic_basis(rep,order,perp)
    matrices=[np.linalg.matrix_power(rep,k)for k in range(order)]
    if dihedral:matrices.extend([r@perp for r in matrices.copy()])
    averaged=sum(r@raw@r.conj().T for r in matrices)/len(matrices)
    basis=layout['basis'];solved=np.zeros_like(raw)
    for k,u in enumerate(basis):
     block=u.conj().T@raw@u
     peer=layout['partner'][k]
     if peer<0:peer=next((p for p,src in enumerate(layout['partner'])if src==k),-1)
     if peer>=0:block=(block+(basis[peer].conj().T@raw@basis[peer]).conj())*.5
     if layout['dtype'][k]==1:block=block.real
     solved+=u@block@u.conj().T
    self.assertLess(np.max(abs(averaged-solved)),1e-12)

 def test_improper_map_schema_declares_axial_lift(self):
  helper=load('generic_export_map_schema',ROOT/'geometry/generic_soc_export.py')
  with tempfile.TemporaryDirectory()as d:
   d=Path(d);inp=d/'soc.dat'
   inp.write_text('''System.Name s
scf.SpinPolarization NC
scf.SpinOrbit.Coupling on
Atoms.Number 1
<Definition.of.Atomic.Species
X X7.0-s1p1 X_PBE19
Definition.of.Atomic.Species>
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
10 0 0
0 10 0
0 0 10
Atoms.UnitVectors>
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.SpeciesAndCoordinates
1 X 0 0 0 1 1 0 0 0 0
Atoms.SpeciesAndCoordinates>
''')
   ops=[]
   for q in [np.eye(3),-np.eye(3)]:
    ops.append({'W':q.tolist(),'t':[0,0,0],'R_cart_orthogonal':q.tolist(),'atom_permutation':[0],'image_shift':[[0,0,0]]})
   result=helper.build(inp,{'operations':ops,'validation':{'geometry_revalidated':True,'symprec_Ang':1e-6}},d/'out.socrs',d/'out.json')
   self.assertEqual(struct.unpack('<8i',(d/'out.socrs').read_bytes()[:32])[:6],(0x52534433,3,1,4,2,1))
   self.assertEqual(result['spin_lift'],'axial_detR_R')
   self.assertLess(result['max_spin_adjoint_group_error'],1e-12)
if __name__=='__main__':unittest.main()
