"""Execute the actual merged DFT one-shot guard for scalar and spinor paths."""
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
class OneShotGuardScope(unittest.TestCase):
 def test_local_hs_does_not_turn_full_nsoc_into_one_step(self):
  source=(ROOT/'overrides/nsoc/source/DFT.unified.c').read_text()
  lines=source.splitlines();at=next(i for i,x in enumerate(lines)if 'CPU_MEMORY_PLAN rejects changing one-step' in x)
  guard=lines[at-1].strip();self.assertTrue(guard.startswith('if(')and guard.endswith('){'))
  expression=guard[3:-2]
  harness='''#include <string.h>\n#include <stdio.h>\nint SpinP_switch,SO_switch,DFTSCF_loop,owner;\nconst char *local_hs;\nint CpuOwnerSparseRequested(void){return owner;}\nint check(void){return EXPR;}\nint main(void){\n SpinP_switch=0;SO_switch=0;DFTSCF_loop=200;owner=0;local_hs="1";if(check())return 11;\n local_hs="0";if(check())return 12;\n SpinP_switch=3;SO_switch=1;DFTSCF_loop=1;owner=1;local_hs="1";if(check())return 13;\n DFTSCF_loop=2;if(!check())return 14;\n owner=0;local_hs="0";if(check())return 15;\n printf("SOC one-shot guard scoped; scalar full-SCF owner storage allowed\\n");return 0;\n}\n'''.replace('EXPR',expression)
  with tempfile.TemporaryDirectory(prefix='soc-guard-') as folder:
   p=Path(folder);(p/'test.c').write_text(harness)
   compile=subprocess.run(['/usr/bin/gcc',str(p/'test.c'),'-o',str(p/'test')],capture_output=True,text=True)
   self.assertEqual(compile.returncode,0,compile.stderr)
   result=subprocess.run([str(p/'test')],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,f'guard returned {result.returncode}; scalar must not inherit SOC one-step restriction')
if __name__=='__main__':unittest.main()
