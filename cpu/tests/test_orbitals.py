"""OpenMX real harmonics: polynomial convention, parity and group laws."""
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

MODULE = Path(__file__).resolve().parents[1] / 'geometry/openmx_orbitals.py'


class TestOrbitals(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not MODULE.exists():
            cls.mod = None
        else:
            spec = importlib.util.spec_from_file_location('local_openmx_orbitals', MODULE)
            cls.mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cls.mod)

    def api(self):
        self.assertIsNotNone(self.mod, 'repository-local orbital provider is missing')
        return self.mod

    def test_comp2real_order_and_unitarity(self):
        tables = self.api().angular_tables()
        self.assertEqual(len(tables), 4)
        for l, c in enumerate(tables):
            self.assertEqual(c.shape, (2*l+1, 2*l+1))
            assert_allclose(c @ c.conj().T, np.eye(len(c)), atol=5e-15)
        # p ordering is px, py, pz, with complex columns m=-1,0,+1.
        assert_allclose(tables[1][2], [0, 1, 0], atol=0)
        tables[0][0, 0] = 4
        self.assertEqual(self.api().angular_tables()[0][0, 0], 1)

    def test_p_is_cartesian_active_rotation(self):
        mod = self.api()
        for seed in range(8):
            r = Rotation.random(random_state=seed).as_matrix()
            for q in (r, -r):
                assert_allclose(mod.get_any_rot_orb_twostep('p', q), q, atol=3e-14)

    def test_group_laws_and_improper_parity(self):
        mod = self.api()
        for shell, l in zip('spdf', range(4)):
            for seed in range(12):
                r, s = Rotation.random(2, random_state=seed).as_matrix()
                if seed % 2:
                    r = -r
                if seed % 3:
                    s = -s
                dr = mod.get_any_rot_orb_twostep(shell, r)
                ds = mod.get_any_rot_orb_twostep(shell, s)
                assert_allclose(dr.T @ dr, np.eye(2*l+1), atol=4e-14)
                assert_allclose(mod.get_any_rot_orb_twostep(shell, r @ s), dr @ ds, atol=4e-14)
            assert_allclose(mod.get_any_rot_orb_twostep(shell, -np.eye(3)), (-1)**l*np.eye(2*l+1), atol=4e-14)

    def test_z_and_arbitrary_c2(self):
        mod = self.api()
        angle = .471
        rz = Rotation.from_rotvec([0, 0, angle]).as_matrix()
        for shell in 'df':
            d = mod.get_any_rot_orb_twostep(shell, rz)
            self.assertAlmostEqual(d[0, 0], 1, places=14)
            for position, m in enumerate(([2, 1] if shell == 'd' else [1, 2, 3])):
                pair = slice(2*position+1, 2*position+3)
                assert_allclose(d[pair, pair], [[np.cos(m*angle), -np.sin(m*angle)], [np.sin(m*angle), np.cos(m*angle)]], atol=3e-14)
        axis = np.array([.2, -.8, .4]); axis /= np.linalg.norm(axis)
        r = 2*np.outer(axis, axis)-np.eye(3)
        for shell in 'spdf':
            d = mod.get_any_rot_orb_twostep(shell, r)
            assert_allclose(d@d, np.eye(len(d)), atol=4e-14)

    def test_explicit_harmonic_polynomials(self):
        mod = self.api()
        rng = np.random.default_rng(142)
        points = rng.normal(size=(80, 3))
        def polynomials(shell, v):
            x, y, z = v.T
            if shell == 's':
                return np.ones((len(v), 1))
            if shell == 'p':
                return v
            if shell == 'd':
                return np.array([(2*z*z-x*x-y*y)/(2*np.sqrt(3)),
                                 (x*x-y*y)/2, x*y, x*z, y*z]).T
            return np.array([z*(2*z*z-3*x*x-3*y*y)/(2*np.sqrt(15)),
                             x*(4*z*z-x*x-y*y)/(2*np.sqrt(10)),
                             y*(4*z*z-x*x-y*y)/(2*np.sqrt(10)),
                             z*(x*x-y*y)/2, x*y*z,
                             x*(x*x-3*y*y)/(2*np.sqrt(6)),
                             y*(3*x*x-y*y)/(2*np.sqrt(6))]).T
        for seed in range(10):
            r = Rotation.random(random_state=seed).as_matrix()
            if seed % 2:
                r = -r
            for shell in 'spdf':
                # f_j(R^-1 x) = sum_i f_i(x) D_ij(R).
                assert_allclose(polynomials(shell, points @ r),
                                polynomials(shell, points) @ mod.get_any_rot_orb_twostep(shell, r),
                                atol=8e-14, rtol=2e-13)

    def test_invalid_and_blocks(self):
        mod = self.api()
        for r in (np.ones((3, 3)), np.eye(2), np.diag([1, 1, np.nan]), np.eye(3)*1.00001):
            with self.assertRaises(ValueError):
                mod.get_any_rot_orb_twostep('p', r)
        with self.assertRaises(ValueError):
            mod.get_any_rot_orb_twostep('g', np.eye(3))
        r = Rotation.from_rotvec([.1, .2, .3]).as_matrix()
        block = mod.orbital_rotation_block('Mo7.0-s3p2d2', r)
        self.assertEqual(block.shape, (19, 19))
        assert_allclose(block[:3, :3], np.eye(3))
        assert_allclose(block[3:6, 3:6], r, atol=3e-14)
        assert_allclose(block.T@block, np.eye(19), atol=4e-14)
        for spec in ('s2p', 's2g1', 's2 s1', 'abc', 'p1s1'):
            with self.assertRaises(ValueError):
                mod.orbital_rotation_block(spec, r)


if __name__ == '__main__':
    unittest.main()
