"""Input-driven symmetry contracts; no native calculation in these tests."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import hashlib

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class InputSymmetry(unittest.TestCase):
    def test_public_boolean_switch_preserves_developer_mode(self):
        m=self.api()
        q=m.request({'scf.symmetry':'on'})
        self.assertEqual(q['mode'],'auto')
        self.assertEqual(q['export_policy'],'symmetrize')
        self.assertEqual(q['hamiltonian_policy'],'symmetrize')
        self.assertEqual(q['tolerance'],5e-5)
        self.assertEqual(q['covariance_tolerance'],1e-5)
        self.assertEqual(m.request({'scf.symmetry':'on','scf.symmetry.mode':'User','scf.symmetry.pointgroup':'C2'})['group'],'C2')
        self.assertEqual(m.request({'scf.symmetry':'off'})['mode'],'off')
        for values in ({'scf.symmetry':'garbage'},{'scf.symmetry':'off','scf.symmetry.mode':'Auto'}):
            with self.assertRaises(ValueError):m.request(values)

    def api(self):
        path = ROOT / 'tools/input_symmetry.py'
        self.assertTrue(path.exists(), 'input-driven symmetry adapter is missing')
        return load('input_symmetry_tests', path)

    def test_auto_request_and_explicit_export(self):
        m = self.api()
        q = m.request({'scf.symmetry.mode': 'Auto', 'scf.symmetry.hsexport': 'Symmetrize'})
        self.assertEqual(q['mode'], 'auto')
        self.assertEqual(q['group'], 'auto')
        self.assertEqual(q['export_policy'], 'symmetrize')

    def test_hamiltonian_projection_is_a_separate_explicit_request(self):
        m = self.api()
        q=m.request({'scf.symmetry.mode':'Auto','scf.symmetry.hsexport':'Symmetrize'})
        self.assertEqual(q.get('hamiltonian_policy'),'free')
        q=m.request({'scf.symmetry.mode':'Auto','scf.symmetry.hsexport':'Symmetrize','scf.symmetry.hamiltonianpolicy':'Symmetrize'})
        self.assertEqual(q['hamiltonian_policy'],'symmetrize')

    def test_user_requires_group_and_valid_tolerance(self):
        m = self.api()
        with self.assertRaisesRegex(ValueError, 'PointGroup'):
            m.request({'scf.symmetry.mode': 'User'})
        for tol in ('nan', '-1', 'inf'):
            with self.assertRaisesRegex(ValueError, 'Tolerance'):
                m.request({'scf.symmetry.mode': 'User', 'scf.symmetry.pointgroup': 'C2',
                           'scf.symmetry.tolerance': tol})

    def test_unimplemented_policy_is_not_ignored(self):
        m = self.api()
        for key in ('geometrypolicy', 'densitypolicy', 'forcepolicy'):
            with self.assertRaises(ValueError):
                m.request({'scf.symmetry.mode': 'Auto', 'scf.symmetry.'+key: 'Project'})
        with self.assertRaisesRegex(ValueError, 'unknown'):
            m.request({'scf.symmetry.mode': 'Auto', 'scf.symmetry.typo': 'on'})

    def test_auto_never_accepts_mismatched_prepared_profiles(self):
        m = self.api()
        q = m.request({'scf.symmetry.mode': 'Auto'})
        with self.assertRaisesRegex(ValueError, 'configure|prepared'):
            m.validate_runtime_request(q, {}, 'nsoc')
        with tempfile.TemporaryDirectory() as d:
            basis = Path(d)/'scalar.nsym'; basis.write_bytes(b'test')
            env = {'OPENMX_CPU_INPUT_SYMMETRY_MODE': 'auto',
                   'OPENMX_CPU_INPUT_SYMMETRY_GROUP': 'D3',
                   'OPENMX_CPU_INPUT_SYMMETRY_TOL': '1e-5',
                   'OPENMX_CPU_INPUT_SYMMETRY_EXPORT': 'raw',
                   'OPENMX_CPU_NSOC_GENERIC_COV_TOL': '1e-10',
                   'OPENMX_CPU_NSOC_SYM': '1', 'OPENMX_CPU_NSOC_BASIS': str(basis)}
            m.validate_runtime_request(q, env, 'nsoc')
            env['OPENMX_CPU_INPUT_SYMMETRY_MODE'] = 'user'
            with self.assertRaises(ValueError):
                m.validate_runtime_request(q, env, 'nsoc')
            env['OPENMX_CPU_INPUT_SYMMETRY_MODE'] = 'auto'
            env['OPENMX_CPU_NSOC_GENERIC_COV_TOL'] = '1e-5'
            with self.assertRaisesRegex(ValueError, 'covariance'):
                m.validate_runtime_request(q, env, 'nsoc')

    def test_off_preserves_legacy_explicit_profiles(self):
        m = self.api()
        m.validate_runtime_request(m.request({}), {'OPENMX_CPU_NSOC_SYM': '1'}, 'nsoc')

    def test_d3_return_optimization_is_not_applied_to_other_groups(self):
        m = self.api()
        self.assertTrue(hasattr(m, 'return_optimizations'), 'missing group-specific capability selection')
        for group in ('C2','D2','Ci','C3','D4'):
            flags = m.return_optimizations(group)
            self.assertEqual(flags['OPENMX_CPU_NSOC_FUSED_RETURN'], '0')
            self.assertEqual(flags['OPENMX_CPU_NSOC_SELECTED_BACK'], '0')
        self.assertEqual(m.return_optimizations('D3')['OPENMX_CPU_NSOC_FUSED_RETURN'], '1')

    def test_declared_export_matches_actual_nsoc_output_mode(self):
        m=self.api()
        with tempfile.TemporaryDirectory()as d:
            basis=Path(d)/'basis';basis.write_bytes(b'test')
            env={'OPENMX_CPU_INPUT_SYMMETRY_MODE':'auto','OPENMX_CPU_INPUT_SYMMETRY_GROUP':'D2',
                 'OPENMX_CPU_INPUT_SYMMETRY_TOL':'1e-5','OPENMX_CPU_INPUT_SYMMETRY_EXPORT':'symmetrize',
                 'OPENMX_CPU_NSOC_GENERIC_COV_TOL':'1e-10','OPENMX_CPU_NSOC_SYM':'1',
                 'OPENMX_CPU_NSOC_BASIS':str(basis),'OPENMX_CPU_REALSPACE_MAP':str(basis),
                 'OPENMX_CPU_REALSPACE_D3':'off'}
            q=m.request({'scf.symmetry.mode':'Auto','scf.symmetry.hsexport':'Symmetrize'})
            with self.assertRaisesRegex(ValueError,'export|output'):
                m.validate_runtime_request(q,env,'nsoc')
            env['OPENMX_CPU_REALSPACE_D3']='export';m.validate_runtime_request(q,env,'nsoc')
            env['OPENMX_CPU_INPUT_SYMMETRY_EXPORT']='raw'
            with self.assertRaisesRegex(ValueError,'export|output'):
                m.validate_runtime_request(m.request({'scf.symmetry.mode':'Auto'}),env,'nsoc')

    def test_generated_artifacts_cannot_be_replaced_after_configure(self):
        m=self.api()
        self.assertTrue(hasattr(m,'artifact_dependencies'),'prepared artifact hashes are not verified')
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'basis';p.write_bytes(b'original')
            config={'symmetry_artifacts_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()}}
            self.assertEqual(m.artifact_dependencies(config),config['symmetry_artifacts_sha256'])
            p.write_bytes(b'replaced')
            with self.assertRaisesRegex(ValueError,'artifact.*changed'):
                m.artifact_dependencies(config)


if __name__ == '__main__':
    unittest.main()
