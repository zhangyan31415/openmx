"""The real native final-log merger retains symmetry near the front of .out."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

CPU = Path(__file__).resolve().parents[1]


class SymmetryOut(unittest.TestCase):
    def test_final_merge_keeps_report_before_input_and_physics_outputs(self):
        source = CPU/'overrides/nsoc/source/Merge_LogFile.symmetry.c'
        self.assertTrue(source.is_file(), 'unified native log merge has no symmetry report hook')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'openmx_common.h').write_text('''
#define YOUSO10 4096
static int mpi_comm_level1=0,Host_ID=0,openmp_threads_num=1,Cnt_switch=0,atomnum=0,SpinP_switch=0,Ngrid1=0,Ngrid2=0;
static char filepath[4096]="./",filename[4096]="sample",Version_OpenMX[]="4.0";
''')
            (root/'mpi.h').write_text('''
static int MPI_Comm_size(int c,int*n){*n=4;return 0;}
static int MPI_Comm_rank(int c,int*r){*r=0;return 0;}
static int MPI_Barrier(int c){return 0;}
''')
            (root/'test.c').write_text('#include "'+str(source)+'"\nint main(void){Merge_LogFile("input.dat");return 0;}\n')
            subprocess.run([shutil.which('cc') or 'cc', '-std=gnu99', '-I'+str(root),
                            str(root/'test.c'), '-o', str(root/'test')],
                           check=True, capture_output=True, text=True)
            (root/'input.dat').write_text('scf.Symmetry on\nAtoms.Number 222\n')
            (root/'sample.DFTSCF').write_text('SCF physical output\n')
            report = root/'report.txt'
            report.write_text('Symmetry analysis -- NSOC\nSelected subgroup     : D3\nGroup operations      : 6\n')
            import os
            env = {**os.environ, 'OPENMX_CPU_INPUT_SYMMETRY_REPORT': str(report)}
            subprocess.run([str(root/'test')], cwd=root, env=env, check=True)
            text = (root/'sample.out').read_text()
            self.assertEqual(text.count('Symmetry analysis -- NSOC'), 1)
            self.assertLess(text.index('Symmetry analysis'), text.index('Atoms.Number'))
            self.assertIn('SCF physical output', text)
            env.pop('OPENMX_CPU_INPUT_SYMMETRY_REPORT')
            subprocess.run([str(root/'test')], cwd=root, env=env, check=True)
            self.assertNotIn('Symmetry analysis', (root/'sample.out').read_text())


if __name__ == '__main__':
    unittest.main()
