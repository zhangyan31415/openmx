"""The distributable CPU target must be self-contained under cpu/."""
import importlib.util
from pathlib import Path
import unittest

CPU = Path(__file__).resolve().parents[1]

def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value)
    return value


class PublicLayout(unittest.TestCase):
    def test_preparation_dependencies_are_inside_cpu(self):
        prep=module('public_prepare',CPU/'geometry/prepare.py')
        self.assertTrue(prep.TOOLS.resolve().is_relative_to(CPU))
        self.assertTrue((CPU/'geometry/prepare_soc_basis.py').is_file())

    def test_metadata_providers_are_repository_local(self):
        for name in ('symmetry_detect','prepare_nsoc_auto','prepare_realspace_d3','neighbor_support'):
            self.assertTrue((CPU/'geometry/metadata'/f'{name}.py').is_file(),name)
        scalar=module('public_d3',CPU/'overrides/nsoc/metadata/prepare_paired_d3.py')
        self.assertTrue(scalar.PARENT_TOOLS.resolve().is_relative_to(CPU))

    def test_public_readme_has_cpu_quickstart(self):
        text=(CPU.parent/'README.md').read_text()
        self.assertIn('cpu/build.py',text)
        self.assertIn('scf.Symmetry on',text)
        self.assertIn('mpirun',text)
        self.assertIn('cpu/README.md',text)
        for private_keyword in ('scf.Symmetry.Mode','scf.Symmetry.PointGroup'):
            self.assertNotIn(private_keyword,text)


if __name__=='__main__':unittest.main()
