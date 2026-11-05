from pathlib import Path
import importlib.util,tempfile,unittest,json,struct
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
TOOLS=ROOT/'cpu/geometry/metadata'
def module(name):
 p=TOOLS/(name+'.py');s=importlib.util.spec_from_file_location('local_'+name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
INPUT='''System.Name prep
DATA.PATH /path/to/DFT_DATA19/
scf.SpinPolarization off
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
12 0 0
0 12 0
0 0 12
Atoms.UnitVectors>
Species.Number 1
<Definition.of.Atomic.Species
H H5.0-s1p1 H_PBE19
Definition.of.Atomic.Species>
Atoms.Number 1
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.SpeciesAndCoordinates
1 H 0 0 0 .5 .5
Atoms.SpeciesAndCoordinates>
'''
class TestPreparation(unittest.TestCase):
 def test_d3_map_accepts_both_hexagonal_c2_axis_families(self):
  d=module('symmetry_detect');r=module('prepare_realspace_d3')
  with tempfile.TemporaryDirectory()as tmp:
   p=Path(tmp)/'input.dat'
   cell=np.array([[12.,0,0],[-6,6*np.sqrt(3),0],[0,0,20.]])
   p.write_text(INPUT.replace('12 0 0\n0 12 0\n0 0 12','\n'.join(' '.join(map(str,row))for row in cell)))
   s=d.parse_openmx(p);theta=2*np.pi/3
   q3=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
   for index,angle in enumerate([0,np.pi/6]):
    axis=np.array([np.cos(angle),np.sin(angle),0]);q2=2*np.outer(axis,axis)-np.eye(3)
    qs=[np.eye(3),q3,q3@q3,q2,q3@q2,q3@q3@q2]
    ws=[np.rint(np.linalg.inv(cell.T)@q@cell.T).astype(int)for q in qs]
    g=d._build_operation_set(s,ws,[np.zeros(3)]*6,1e-6,1.25,source='test_rotated_d3')
    self.assertTrue(g['valid'])
    out=r.generate(p,Path(tmp)/f'd3_{index}.rsd3',operations_metadata={'selected':g},atom_tol=1e-5)
    self.assertEqual(out['operation_count'],6)

 def test_orbital_provider_is_repository_local(self):
  d=module('symmetry_detect');f=d._load_orbital_rotation_provider()
  self.assertIsNotNone(f)
  self.assertTrue(Path(f.__code__.co_filename).resolve().is_relative_to(ROOT/'cpu/geometry'))
 def test_unwrapped_mapping_records_original_image(self):
  d=module('symmetry_detect')
  with tempfile.TemporaryDirectory()as tmp:
   p=Path(tmp)/'input.dat';p.write_text(INPUT.replace('1 H 0 0 0','1 H 120 -96 0'))
   s=d.parse_openmx(p);op=d._map_operation(s,np.diag([-1,-1,1]),np.zeros(3),1e-8)
   self.assertEqual(op['perm'].tolist(),[0]);self.assertEqual(op['image_shift'].tolist(),[[-20,16,0]])
 def test_realspace_preserves_validated_subgroup(self):
  d=module('symmetry_detect');r=module('prepare_realspace_d3')
  with tempfile.TemporaryDirectory()as tmp:
   p=Path(tmp)/'input.dat';p.write_text(INPUT);s=d.parse_openmx(p)
   g=d._build_operation_set(s,[np.eye(3,dtype=int),np.diag([-1,-1,1])],[np.zeros(3),np.zeros(3)],1e-6,1.25,source='unit')
   self.assertTrue(g['valid'])
   out=Path(tmp)/'c2.rsd3'
   info=r.generate(p,out,symprec=1e-6,atom_tol=1e-5,include_improper=True,operations_metadata={'selected':g})
   self.assertEqual(info['operation_count'],2)
   with out.open('rb')as f:self.assertEqual(struct.unpack('<8i',f.read(32))[4],2)
   self.assertTrue(info['geometry_modified']is False)
 def test_near_hexagonal_group_uses_one_invariant_metric(self):
  d=module('symmetry_detect')
  with tempfile.TemporaryDirectory()as tmp:
   text=INPUT.replace('12 0 0\n0 12 0\n0 0 12','12.0001 0 0\n-6 10.392304845413264 0\n0 0 20')
   p=Path(tmp)/'input.dat';p.write_text(text);s=d.parse_openmx(p)
   w=np.array([[0,-1,0],[1,-1,0],[0,0,1]])
   g=d._build_operation_set(s,[np.eye(3,dtype=int),w,w@w],[np.zeros(3)]*3,.01,1.25,source='unit')
   self.assertTrue(g['valid'])
   self.assertLess(g['orbital_representation']['max_composition_error'],2e-13)
if __name__=='__main__':unittest.main()
