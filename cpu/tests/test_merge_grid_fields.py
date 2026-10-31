import json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parent
class UnifiedGridFieldHeaders(unittest.TestCase):
    def test_all_selected_consumers_use_one_header(self):
        soc=json.loads((ROOT/'overrides/soc/SOURCES.json').read_text())
        selection=json.loads((ROOT/'overrides/nsoc/UNIFIED_SELECTION.json').read_text())
        consumers=[]
        for unit in soc['units']+selection['extra_units']:
            source=REPO/selection['replacements'].get(unit['object'],unit['source'])
            text=source.read_text()
            self.assertNotIn('#include "cpu_grid_fields_shared.h"',text,str(source))
            if '#include "cpu_grid_fields_shared.unified.h"' in text:consumers.append(unit['object'])
        self.assertEqual(len(consumers),7)
    def test_layout_is_agreed_before_grid_allocation(self):
        source=(ROOT/'overrides/nsoc/source/truncation.unified.c').read_text()
        self.assertLess(source.index('CpuGridFields_Configure('),source.index('CpuGridFields_Init('))
        self.assertIn('SpinP_switch==3?4:2',source)
        header=(ROOT/'overrides/nsoc/source/cpu_grid_fields_shared.unified.h').read_text()
        self.assertIn('#define CPU_GF_FIELDS (3*cpu_grid_fields_width+2)',header)
        self.assertIn('generation[CPU_GF_MAX_FIELDS]',(ROOT/'overrides/nsoc/source/cpu_grid_fields_shared.unified.c').read_text())
if __name__=='__main__':unittest.main()
