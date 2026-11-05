"""End-to-end CLI includes NC geometry metadata and scalar AO revalidation."""
from pathlib import Path
import unittest,subprocess,tempfile,sys,json,hashlib,os
SCRIPT=Path(__file__).resolve().parents[1]/'geometry/prepare_soc_basis.py'
class AutoCliTests(unittest.TestCase):
    def fixture(self,root):
        inp=root/'nc.dat';inp.write_text('''Species.Number 1
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
scf.SpinPolarization NC
scf.SpinOrbit.Coupling on
''');return inp
    def invoke(self,args):
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
        return subprocess.run([sys.executable,str(SCRIPT),*map(str,args)],capture_output=True,text=True,env=env)
    def test_auto_nc_input_and_explicit_json_compatibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inp=self.fixture(root);before=hashlib.sha256(inp.read_bytes()).hexdigest();out=root/'auto.gsoc'
            run=self.invoke(['--input',inp,'--output',out,'--geometry','strict','--symprec-ladder','1e-6','3e-6','1e-5'])
            self.assertEqual(run.returncode,0,run.stderr)
            self.assertTrue(out.exists());meta=json.loads((root/'auto.symmetry.json').read_text());info=json.loads((root/'auto.json').read_text())
            self.assertEqual(meta['spin_mode'],'spinor');self.assertFalse(meta['orbital_representation_checked'])
            self.assertTrue(info['validation']['geometry_revalidated']);self.assertTrue(info['validation']['scalar_ao_representation_revalidated']);self.assertEqual(info['selected_subgroup'],'D2')
            self.assertEqual(hashlib.sha256(inp.read_bytes()).hexdigest(),before)
            explicit=self.invoke(['--input',inp,'--output',root/'explicit.gsoc','--symmetry-json',root/'auto.symmetry.json'])
            self.assertEqual(explicit.returncode,0,explicit.stderr)
            self.assertEqual(out.read_bytes(),(root/'explicit.gsoc').read_bytes())
    def test_refuse_overwrite_and_bad_ladder(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inp=self.fixture(root);out=root/'protected.gsoc';out.write_bytes(b'preserved reference')
            run=self.invoke(['--input',inp,'--output',out])
            self.assertNotEqual(run.returncode,0);self.assertIn('already exists',run.stderr)
            self.assertEqual(out.read_bytes(),b'preserved reference')
            run=self.invoke(['--input',inp,'--output',root/'bad.gsoc','--symprec-ladder','1e-5','1e-6'])
            self.assertNotEqual(run.returncode,0);self.assertIn('increasing',run.stderr)
if __name__=='__main__':unittest.main()
