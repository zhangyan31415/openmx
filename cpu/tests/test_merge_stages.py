import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('merge_stages', ROOT/'overrides/merge_stages.py')

class MergeStagesTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SPEC.origin and pathlib.Path(SPEC.origin).exists(), 'context-checked stage merger missing')
        self.module = importlib.util.module_from_spec(SPEC)
        SPEC.loader.exec_module(self.module)

    def test_preserves_independent_soc_insertions(self):
        base = 'header\na\nb\nc\nd\ne\nf\ng\n'
        nsoc = 'header\na\nb\nNSOC\nc\nd\ne2\nf\ng\n'
        soc = 'header\nSOC guard\na\nb\nc\nd\ne\nf\ng\nSOC release\n'
        actual = self.module.merge_text(base, nsoc, soc)
        self.assertEqual(actual, 'header\nSOC guard\na\nb\nNSOC\nc\nd\ne2\nf\ng\nSOC release\n')

    def test_already_shared_delta_is_not_duplicated(self):
        self.assertEqual(self.module.merge_text('a\nb\nc\n', 'a\nNSOC\nb\nc\n', 'SOC\na\nNSOC\nb\nc\n'), 'SOC\na\nNSOC\nb\nc\n')

    def test_rejects_changed_source_context(self):
        with self.assertRaises(ValueError):
            self.module.merge_text('a\nb\nc\n', 'a\nx\nc\n', 'a\nchanged\nc\n')

    def test_rejects_ambiguous_insertion(self):
        with self.assertRaises(ValueError):
            self.module.merge_text('a\nb\n', 'a\nNSOC\nb\n', 'a\nb\na\nb\n')

    def test_unified_truncation_keeps_nsoc_and_soc_memory_contracts(self):
        nsoc = (ROOT/'overrides/nsoc/source/truncation.symhistory.c').read_text()
        soc = (ROOT/'overrides/soc/source/truncation.generic-soc.c').read_text()
        actual = self.module.merge_truncation(nsoc, soc)
        for marker in ['CpuNSOCGhostAgree', 'CpuLocalHSPlan', 'CpuLocalOLPPlan', 'CpuOneStepHistoryPlan',
                       'CpuNSocSymHistoryEnabled', 'CpuNSOCImagSlots', 'CpuNSocHamiltonianHistorySlots']:
            self.assertIn(marker, actual)
        self.assertEqual(actual.count('int cpu_local_hs='), 1)
        self.assertEqual(actual.count('int cpu_local_olp='), 0)

    def test_real_dft_merge_retains_both_lifecycles(self):
        repo = ROOT.parent
        stock = repo/'cpu/vendor/openmx4.0/source/DFT.c'
        if not stock.exists():
            stock = pathlib.Path('snapshots/home/software/openmx4.0/source/DFT.c')
        nsoc = ROOT/'overrides/nsoc/source/DFT.direct.c'
        soc = repo/'local-records/cpu_d3/generic_soc/integration/DFT.compact.c'
        self.assertTrue(nsoc.exists(), 'recovered NSOC DFT source missing')
        merged = self.module.merge_text(stock.read_text(), nsoc.read_text(), soc.read_text())
        scalar_start=merged.index('void Allocate_Free_Cluster_Col(int todo_flag)\n{')
        scalar_end=merged.index('void Allocate_Free_Cluster_NonCol(int todo_flag)\n{',scalar_start)
        self.assertTrue(scalar_start < merged.index('CpuNSym_DirectWorkspaceRecord') < scalar_end, 'NSOC hook escaped scalar allocator')
        for marker in ['CpuNSym_ReserveEigenScratch', 'CpuNSym_ConvergenceBlocked',
                       'CpuGSOC_ValidateRuntime', 'CpuGSOC_CompactRelease',
                       'CPU_NSOC_SCF_BUFFERS_RELEASED', 'CpuC3_ReleaseDenseOneShot']:
            self.assertIn(marker, merged)

if __name__ == '__main__':
    unittest.main()
