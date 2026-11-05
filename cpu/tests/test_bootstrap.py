"""Automatic native-start metadata preparation; never launch OpenMX."""
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

CPU = Path(__file__).resolve().parents[1]


def helper():
    path = CPU/'geometry/bootstrap.py'
    if not path.exists():
        raise AssertionError('native automatic preparation helper is missing')
    spec = importlib.util.spec_from_file_location('test_native_bootstrap', path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


INPUT = '''System.Name native
DATA.PATH data
scf.Symmetry on
scf.Symmetry.Mode User
scf.Symmetry.PointGroup D2
scf.EigenvalueSolver Cluster
scf.eigen.lib elpa2
scf.SpinPolarization off
scf.XcType GGA-PBE
scf.energycutoff 150
scf.ElectronicTemperature 300
scf.Mixing.Type RMM-DIIS
scf.Mixing.History 30
scf.Mixing.StartPulay 15
scf.maxIter 40
HS.fileout off
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
12 0 0
0 14 0
0 0 16
Atoms.UnitVectors>
Species.Number 1
<Definition.of.Atomic.Species
H H2.0-s1p1 H_PBE19
Definition.of.Atomic.Species>
Atoms.Number 1
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.SpeciesAndCoordinates
1 H 0 0 0 .5 .5
Atoms.SpeciesAndCoordinates>
'''


class Bootstrap(unittest.TestCase):
    def test_simple_mixing_does_not_activate_incompatible_memory_paths(self):
        m=helper()
        for text in (INPUT.replace('RMM-DIIS','Simple'), INPUT.replace('scf.Mixing.Type RMM-DIIS\n','')):
            with tempfile.TemporaryDirectory() as d:
                path,out=self.fixture(Path(d),text)
                result=m.bootstrap(path,out,'nsoc',out/'environment.txt',runtime_cwd=Path(d))
                for key in ('LOCAL_HS','LOCAL_OLP','GRID_FIELDS_SHARED','VNA_SHARED','NSOC_RELEASE_SCF_BUFFERS','NSOC_OWNER_DM'):
                    self.assertNotEqual(result['environment'].get('OPENMX_CPU_'+key), '1', key)

    def test_symmetry_history_preserves_non30_user_history(self):
        m=helper();api=m.module('history_api',CPU/'tools/input_symmetry.py');q=api.request({'scf.symmetry':'on'})
        prepared={'selected_group':'D3','SOC_route':'legacy','NSOC_environment_additions':{'OPENMX_CPU_REALSPACE_D3':'scf'},'SOC_environment_additions':{}}
        for history in ('15','30','35',None):
            values={'scf.eigenvaluesolver':'Cluster','scf.mixing.type':'rmm-diis'}
            if history is not None:values['scf.mixing.history']=history
            env=m.profiles(prepared,q,api,values=values)['nsoc']
            self.assertEqual(env.get('OPENMX_CPU_NSOC_SYM_HISTORY','0'),'1' if history=='30' else '0')

    def test_native_zeeman_key_keeps_specialized_soc_disabled(self):
        m=helper()
        self.assertEqual(m.preferred_soc_route({'scf.eigenvaluesolver':'cluster','scf.mixing.type':'rmm-diis','scf.nc.zeeman.spin':'on'}),'generic')
    def test_memory_optimizations_are_automatic_for_scalar_cluster(self):
        m = helper()
        with tempfile.TemporaryDirectory() as d:
            path, out = self.fixture(Path(d))
            env = m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=Path(d))['environment']
            for name in ('NSOC_ASSEMBLY_WORKSPACE', 'NSOC_OWNER_DM', 'NSOC_IMAG_COMPACT',
                         'NSOC_COMPACT_WORK1', 'NSOC_H_HISTORY_COMPACT', 'NSOC_EIGEN_SCRATCH',
                         'GRID_SHARED', 'GRID_FIELDS_SHARED', 'LOCAL_HS', 'LOCAL_OLP',
                         'NSOC_RELEASE_SCF_BUFFERS', 'H_GRID_BLAS3'):
                self.assertEqual(env['OPENMX_CPU_'+name], '1', name)
            self.assertEqual(env['OPENMX_CPU_NSOC_DIRECT_SECTORS'], '0')

    def test_threaded_profile_preserves_memory_savings_without_single_thread_caches(self):
        m = helper()
        with tempfile.TemporaryDirectory() as d:
            path, out = self.fixture(Path(d))
            result = m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=Path(d), threads=4)
            env = result['environment']
            self.assertEqual(env['OPENMX_CPU_GRID_SHARED'], '1')
            self.assertEqual(env['OPENMX_CPU_NSOC_OWNER_DM'], '1')
            self.assertEqual(env['OPENMX_CPU_GRID_NUMA'], '0')
            self.assertEqual(env['OPENMX_CPU_VNA_SHARED'], '0')
            self.assertEqual(result['parent_threads'], 4)

    def test_numa_replica_cache_requires_uniform_single_cpu_rank_binding(self):
        m=helper()
        for pinned in (False, True):
            with self.subTest(pinned=pinned), tempfile.TemporaryDirectory() as d:
                path,out=self.fixture(Path(d))
                result=m.bootstrap(path,out,'nsoc',out/'environment.txt',runtime_cwd=Path(d),ranks_single_cpu_pinned=pinned)
                self.assertEqual(result['environment']['OPENMX_CPU_GRID_NUMA'], '1' if pinned else '0')

    def test_auto_soc_route_selects_verified_hexagonal_path_only(self):
        m = helper()
        import numpy as np
        values = {'scf.eigenvaluesolver': 'Cluster', 'scf.mixing.type':'rmm-diis'}
        cell = np.array([[12., 0, 0], [-6, 6*np.sqrt(3), 0], [0, 0, 20.]])
        self.assertEqual(m.preferred_soc_route(values), 'auto')
        self.assertEqual(m.preferred_soc_route({**values, 'scf.symmetry.covariancetolerance': '1e-6'}), 'generic')
        prep=m.module('route_prepare', CPU/'geometry/prepare.py')
        structure={'cell':cell}
        self.assertEqual(prep.automatic_soc_route('D3', structure, 1e-5), 'legacy')
        self.assertEqual(prep.automatic_soc_route('D2', structure, 1e-5), 'generic')
        self.assertEqual(prep.automatic_soc_route('D3', {'cell':np.eye(3)}, 1e-5), 'generic')

    def test_native_d3_defaults_use_complete_legacy_soc_profile(self):
        m = helper()
        text = INPUT.replace('PointGroup D2', 'PointGroup D3').replace('0 14 0', '-6 10.392304845413264 0')
        text = text.replace('scf.SpinPolarization off', 'scf.SpinPolarization nc\nscf.SpinOrbit.Coupling on\nscf.restart c2n').replace('scf.maxIter 40', 'scf.maxIter 1')
        with tempfile.TemporaryDirectory() as d:
            path, out = self.fixture(Path(d), text)
            result = m.bootstrap(path, out, 'soc', out/'environment.txt', runtime_cwd=Path(d))
            env = result['environment']
            self.assertEqual(result['SOC_route'], 'legacy')
            for key in ('GRID_SHARED','GRID_FIELDS_SHARED','SOC_COMPACT','OWNER_SPARSE','SOC_TR','SOC_REAL','HNL_BATCH'):
                self.assertEqual(env['OPENMX_CPU_'+key], '1', key)
            self.assertEqual(env['OPENMX_CPU_C3'], '3')
            self.assertEqual(env['OPENMX_CPU_GSOC'], '0')

    def test_native_logical_aliases_match_stock_inputtools(self):
        m = helper()
        for value in ('on', 'yes', 'true', '.true.', 'ok', 'ON'):
            with self.subTest(value=value):
                self.assertTrue(m.native_logical(value))
                m.check_scope({'scf.eigenvaluesolver': 'Cluster', 'scf.eigen.lib': 'elpa2',
                               'scf.spinpolarization': 'nc', 'scf.spinorbit.coupling': value,
                               'scf.restart': 'c2n', 'scf.maxiter': '1'}, {}, 'soc')
        for value in ('off', 'no', 'false', '.false.', 'ng', 'OFF'):
            with self.subTest(value=value):
                self.assertFalse(m.native_logical(value))
        for value in ('1', '0', 'invalid'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                m.native_logical(value)

    def test_soc_native_logical_export_alias(self):
        m = helper()
        text = INPUT.replace('scf.SpinPolarization off', 'scf.SpinPolarization nc\nscf.SpinOrbit.Coupling .true.\nscf.restart c2n').replace('scf.maxIter 40', 'scf.maxIter 1').replace('HS.fileout off', 'HS.fileout ok')
        with tempfile.TemporaryDirectory() as d:
            path, out = self.fixture(Path(d), text)
            result = m.bootstrap(path, out, 'soc', out/'environment.txt', runtime_cwd=Path(d))
            self.assertTrue(result['actual_HS_export_requested'])
            self.assertEqual(result['environment']['OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE'], '1')

    def test_soc_without_hs_output_projects_solve_without_enabling_export_hook(self):
        m=helper()
        text=INPUT.replace('scf.SpinPolarization off','scf.SpinPolarization nc\nscf.SpinOrbit.Coupling on\nscf.restart c2n').replace('scf.maxIter 40','scf.maxIter 1')
        with tempfile.TemporaryDirectory()as d:
            path,out=self.fixture(Path(d),text)
            result=m.bootstrap(path,out,'soc',out/'environment.txt',runtime_cwd=Path(d))
            self.assertEqual(result['environment']['OPENMX_CPU_GSOC_PROJECT'],'1')
            self.assertEqual(result['environment']['OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE'],'0')
    def fixture(self, root, text=INPUT):
        (root/'data/PAO').mkdir(parents=True)
        (root/'data/PAO/H2.0.pao').write_text('radial.cutoff.pao 2.0\n')
        (root/'inputs').mkdir()
        path = root/'inputs/input.dat'
        path.write_text(text)
        out = root/'private'
        out.mkdir()
        return path, out

    def test_single_input_generates_checked_d2_profile_without_rewriting_input(self):
        m = helper()
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, out = self.fixture(root)
            result = m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=root)
            self.assertEqual(path.read_text(), INPUT)
            self.assertEqual(result['selected_group'], 'D2')
            env = dict(line.split('=', 1) for line in (out/'environment.txt').read_text().splitlines())
            self.assertEqual(env['OPENMX_CPU_INPUT_SYMMETRY_MODE'], 'user')
            self.assertEqual(env['OPENMX_CPU_INPUT_SYMMETRY_GROUP'], 'D2')
            self.assertEqual(env['OPENMX_CPU_REALSPACE_D3'], 'scf')
            self.assertEqual(env['OPENMX_CPU_NSOC_FUSED_RETURN'], '0')
            self.assertEqual(env['OPENMX_CPU_NSOC_SELECTED_BACK'], '0')
            self.assertTrue(all(k.startswith('OPENMX_CPU_') for k in env))
            self.assertEqual(float(env['OPENMX_CPU_NSOC_GENERIC_COV_TOL']), 1e-5)
            for stage in ('nsoc', 'soc'):
                text = (out/'canonical'/f'{stage}.dat').read_text()
                self.assertIn('METADATA ONLY', text)
                self.assertIn(str(root/'data'), text)
                self.assertIn('scf.Mixing.History 30', text)
                self.assertIn('scf.Mixing.StartPulay 15', text)
            provenance = json.loads((out/'BOOTSTRAP.json').read_text())
            self.assertEqual(provenance['actual_input'], str(path))
            self.assertFalse(provenance['actual_input_modified'])
            self.assertFalse(provenance['native_calculation_started'])

    def test_soc_uses_same_group_and_symmetrized_generic_export(self):
        m = helper()
        text = INPUT.replace('scf.SpinPolarization off', 'scf.SpinPolarization NC\nscf.SpinOrbit.Coupling on\nscf.Restart C2N').replace('scf.maxIter 40', 'scf.maxIter 1').replace('HS.fileout off','HS.fileout on')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, out = self.fixture(root, text)
            result = m.bootstrap(path, out, 'soc', out/'environment.txt', runtime_cwd=root)
            self.assertEqual(result['selected_group'], 'D2')
            env = dict(line.split('=', 1) for line in (out/'environment.txt').read_text().splitlines())
            self.assertEqual(env['OPENMX_CPU_GSOC'], '1')
            self.assertEqual(env['OPENMX_CPU_GSOC_PROJECT'], '1')
            self.assertEqual(env['OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE'], '1')
            self.assertTrue(Path(env['OPENMX_CPU_SOC_NATIVE_MAP']).is_file())

    def test_simple_on_selects_automatically_and_logs_only_summary(self):
        m = helper()
        text = INPUT.replace('scf.Symmetry.Mode User\nscf.Symmetry.PointGroup D2\n', '')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, out = self.fixture(root, text)
            captured = io.StringIO()
            with contextlib.redirect_stdout(captured):
                result = m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=root)
            self.assertEqual(result['request']['mode'], 'auto')
            self.assertEqual(result['selected_group'], 'D2')
            self.assertLess(len(captured.getvalue()), 1024)
            self.assertEqual(len(captured.getvalue().splitlines()), 1)
            self.assertTrue(Path(result['preparation_log']).is_file())

    def test_failed_metadata_preparation_leaves_no_success_environment(self):
        m = helper()
        text = INPUT.replace('scf.Symmetry.PointGroup D2', 'scf.Symmetry.PointGroup D3')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, out = self.fixture(root, text)
            with self.assertRaises(ValueError):
                m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=root)
            self.assertFalse((out/'environment.txt').exists())
            self.assertEqual(json.loads((out/'BOOTSTRAP.json').read_text())['status'], 'failed')
            self.assertEqual(path.read_text(), text)

    def test_unsupported_kpoints_fail_without_final_environment_file(self):
        m = helper()
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, out = self.fixture(root, INPUT+'scf.Kgrid 2 1 1\n')
            with self.assertRaisesRegex(ValueError, 'Gamma'):
                m.bootstrap(path, out, 'nsoc', out/'environment.txt', runtime_cwd=root)
            self.assertFalse((out/'environment.txt').exists())

    def test_rejects_multiple_soc_steps_and_magnetic_populations(self):
        m = helper()
        for change, stage, pattern in [
            (INPUT.replace('scf.SpinPolarization off', 'scf.SpinPolarization NC\nscf.SpinOrbit.Coupling on\nscf.Restart C2N'), 'soc', 'one.step'),
            (INPUT.replace('.5 .5', '.6 .4'), 'nsoc', 'nonmagnetic'),
        ]:
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                path, out = self.fixture(root, change)
                with self.assertRaisesRegex(ValueError, pattern):
                    m.bootstrap(path, out, stage, out/'environment.txt', runtime_cwd=root)
                self.assertFalse((out/'environment.txt').exists())


if __name__ == '__main__':
    unittest.main()
