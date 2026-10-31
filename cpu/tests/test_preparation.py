from pathlib import Path
import importlib.util,tempfile,unittest,json,struct
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
TOOLS=ROOT/'local-records/cpu_d3/historical/geometry_tools'
def module(name):
 p=TOOLS/(name+'.py');s=importlib.util.spec_from_file_location('local_'+name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
INPUT='''System.Name prep
DATA.PATH /data/home/apps/openmx/openmx3.9/DFT_DATA19/
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
