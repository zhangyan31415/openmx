import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
NSOC=ROOT/'overrides/nsoc'

class UnifiedMemoryCoverage(unittest.TestCase):
    def test_positive_history_header_preserves_original_archive(self):
        original=(NSOC/'source/cpu_nsoc_history_capacity.h').read_text()
        self.assertIn('Num_Mixing_pDM==30',original)
        newer=NSOC/'source/cpu_nsoc_history_capacity.unified.h'
        self.assertTrue(newer.exists(),'later verified positive-history header missing')
        self.assertIn('Num_Mixing_pDM>0',newer.read_text())
        for name in ['truncation','Free_Arrays','Mixing_H']:
            text=(NSOC/'source'/f'{name}.unified.c').read_text()
            self.assertIn('#include "cpu_nsoc_history_capacity.unified.h"',text)
            self.assertNotIn('#include "cpu_nsoc_history_capacity.h"',text)

    def test_xc_source_has_exact_original_patch_output(self):
        script=NSOC/'prepare_xc_workspace.py'
        self.assertTrue(script.exists(),'verified XC workspace patch missing')
        spec=importlib.util.spec_from_file_location('xc_patch',script)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        stock=(ROOT/'vendor/openmx4.0/source/Set_XC_Grid.c').read_text()
        self.assertEqual(module.prepare(stock),(NSOC/'source/Set_XC_Grid.workspace.c').read_text())
        header=(NSOC/'source/cpu_xc_workspace.h').read_text()
        self.assertIn('!s||!strcmp(s,"0")?0',header)
        self.assertIn('SpinP_switch!=0||SO_switch!=0||MD_switch!=0',header)

    def test_selection_declares_memory_capabilities(self):
        selection=json.loads((NSOC/'UNIFIED_SELECTION.json').read_text())
        self.assertTrue(selection['capabilities']['nsoc_symmetry'])
        self.assertTrue(selection['capabilities']['nsoc_positive_history'])
        self.assertTrue(selection['capabilities']['nsoc_xc_workspace'])
        self.assertEqual(selection['replacements']['Set_XC_Grid.o'],'cpu/overrides/nsoc/source/Set_XC_Grid.workspace.c')

if __name__=='__main__':unittest.main()
