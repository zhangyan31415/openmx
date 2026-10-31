import importlib.util
import pathlib
import unittest
import numpy as np

ROOT=pathlib.Path(__file__).resolve().parents[1]


def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'overrides/nsoc/metadata'/f'{name}.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LocalPairedGeneratorsTest(unittest.TestCase):
    def test_c3_module_imports_without_external_pythonpath(self):
        module=load('prepare_paired_c3')
        self.assertTrue(module.DEFAULT_HELPER.is_file())
        self.assertNotIn('/results/',str(module.DEFAULT_HELPER))

    def test_d3_basis_obeys_scalar_rotation(self):
        module=load('prepare_paired_d3')
        rotation=np.array([[0,0,1],[1,0,0],[0,1,0]],float)
        c2=np.array([[1,0,0],[0,0,1],[0,1,0]],float)
        basis,sizes,error=module.local_basis(rotation,c2)
        np.testing.assert_allclose(basis.T@basis,np.eye(3),atol=1e-13)
        self.assertEqual(sum(sizes),3)
        self.assertLess(max(error.values()),1e-11)

if __name__=='__main__':unittest.main()
