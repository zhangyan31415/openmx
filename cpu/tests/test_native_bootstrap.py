"""Execute the native bootstrap and adapter without a native SCF calculation.

Only MPI collectives and the large original parser are replaced. The checked-in
adapter, bootstrap C implementation, upstream keyword parser and child process
launch all execute in these tests.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'overrides/nsoc/source'
UPSTREAM = ROOT / 'vendor/openmx4.0/source'


STUB = r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#define YOUSO10 4096
#define MPI_INT 1
#define MPI_CHAR 2
#define MPI_MAX 3
#define MPI_MIN 4
static int mpi_comm_level1=0;
static int openmp_threads_num=1;
static int stress_flag=0,StressConsistency_flag=0;
static int MPI_Comm_rank(int c,int *r){*r=0;return 0;}
static int MPI_Allreduce(const void *a,void *b,int n,int t,int op,int c){memcpy(b,a,n*sizeof(int));return 0;}
static int MPI_Bcast(void *a,int n,int t,int r,int c){return 0;}
static int MPI_Abort(int c,int status){exit(status);return 0;}
#include "Inputtools.h"
static const char *value(const char *key){const char *v=getenv(key);return v?v:"unset";}
void Input_std(char *file){
 printf("STOCK_PARSER NSOC_SYM=%s GSOC=%s MODE=%s TEST=%s PROJECT=%s EXPORT=%s\n",
  value("OPENMX_CPU_NSOC_SYM"),value("OPENMX_CPU_GSOC"),
  value("OPENMX_CPU_INPUT_SYMMETRY_MODE"),value("OPENMX_CPU_TEST"),
  value("OPENMX_CPU_GSOC_PROJECT"),value("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE"));
}
'''

PROTOCOL = r'''
#include <assert.h>
#include "cpu_input_bootstrap.h"
int main(void){
 char error[512];
 assert(cpu_bootstrap_environment("OPENMX_CPU_TEST=/a path=value\n",0,error,sizeof(error))==0);
 assert(cpu_bootstrap_environment("OPENMX_CPU_TEST=/a path=value\n",1,error,sizeof(error))==0);
 assert(!strcmp(getenv("OPENMX_CPU_TEST"),"/a path=value"));
 assert(cpu_bootstrap_environment("PATH=/tmp\n",0,error,sizeof(error))!=0);
 assert(cpu_bootstrap_environment("OPENMX_CPU_A=1\nOPENMX_CPU_A=2\n",0,error,sizeof(error))!=0);
 assert(cpu_bootstrap_environment("OPENMX_CPU_A=bad\rvalue\n",0,error,sizeof(error))!=0);
 assert(cpu_bootstrap_environment("OPENMX_CPU_A=bad\tvalue\n",0,error,sizeof(error))!=0);
 assert(cpu_bootstrap_environment("OPENMX_CPU_BAD-KEY=value\n",0,error,sizeof(error))!=0);
 assert(cpu_bootstrap_environment("",0,error,sizeof(error))!=0);
 setenv("OPENMX_CPU_STALE","1",1);
 setenv("OPENMX_CPU_SYMMETRY_HELPER","helper with spaces",1);
 setenv("OPENMX_CPU_SYMMETRY_PYTHON","python",1);
 assert(cpu_bootstrap_clear_profile(error,sizeof(error))==0);
 assert(!getenv("OPENMX_CPU_STALE"));
 assert(!strcmp(getenv("OPENMX_CPU_SYMMETRY_HELPER"),"helper with spaces"));
 assert(!strcmp(getenv("OPENMX_CPU_SYMMETRY_PYTHON"),"python"));
 puts("PROTOCOL_PASS");return 0;
}
'''

HELPER = r'''
import argparse
from pathlib import Path
p=argparse.ArgumentParser()
for key in ('input','output','stage','environment-file'):
 p.add_argument('--'+key,required=True)
p.add_argument('--threads',type=int,default=1)
p.add_argument('--rank-pinned',type=int,choices=(0,1),default=0)
a=p.parse_args()
print('PARENT_THREADS='+str(a.threads),flush=True)
assert Path(a.input).is_absolute()
values={}
for line in Path(a.input).read_text().splitlines():
 words=line.split()
 if len(words)>1:values[words[0].lower()]=words[1]
print('BOOTSTRAP_CHILD_ONCE',flush=True)
print('CPU_SYMMETRY_SUMMARY stage='+a.stage+' selected_group=C2 operations=2',flush=True)
print('CPU_SYMMETRY_OPERATION index=1/2 W=[1,0,0;0,1,0;0,0,1] t_frac=[0,0,0]',flush=True)
if values.get('test.helper')=='fail':
 raise SystemExit(7)
if values.get('test.helper')=='oversize':
 Path(a.environment_file).write_bytes(b'x'*(16*1024*1024+1))
 raise SystemExit(0)
basis=Path(a.output)/'basis with spaces'
basis.write_text('isolated parser fixture')
env={
 'INPUT_SYMMETRY_MODE':values.get('scf.symmetry.mode','auto').lower(),
 'INPUT_SYMMETRY_GROUP':values.get('scf.symmetry.pointgroup','C2'),
 'INPUT_SYMMETRY_TOL':'5e-5',
 'INPUT_SYMMETRY_EXPORT':'symmetrize',
 'TEST':'/a path=value',
}
if a.stage=='nsoc':
 env.update(NSOC_GENERIC_COV_TOL='1e-5',NSOC_SYM='1',NSOC_BASIS=str(basis),
            REALSPACE_MAP=str(basis),REALSPACE_D3='scf')
else:
 assert a.stage=='soc'
 output_enabled=values.get('hs.fileout','off').lower() in ('on','true','.true.','yes','ok')
 env.update(GSOC_COV_TOL='1e-5',GSOC='1',GSOC_BASIS=str(basis),
            GSOC_EXPORT_SYMMETRIZE='1' if output_enabled else '0',GSOC_PROJECT='1')
 if values.get('test.route')=='legacy':
  env.update(GSOC='0',C3='3',D3='1',INPUT_SYMMETRY_GROUP='D3',
             INPUT_SYMMETRY_COV_TOL='1e-5',C3_BASIS=str(basis),
             SOC_NATIVE_MAP=str(basis),SOC_NATIVE_GROUP='D3',SOC_TR='1',SOC_REAL='1')
Path(a.environment_file).write_text(''.join('OPENMX_CPU_'+key+'='+v+'\n' for key,v in env.items()))
'''


class NativeBootstrap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('cc') or shutil.which('gcc')
        if not compiler:
            raise unittest.SkipTest('a host C compiler is needed for native adapter tests')
        cls.directory = tempfile.TemporaryDirectory(prefix='openmx native bootstrap ')
        cls.build = Path(cls.directory.name)
        (cls.build/'stock_stub.h').write_text(STUB)
        adapter = (SOURCE/'Input_std.symmetry.c').read_text()
        original = '#include "../../../vendor/openmx4.0/source/Input_std.c"'
        if adapter.count(original) != 1:
            raise AssertionError('native adapter upstream include is not uniquely identifiable')
        adapter = adapter.replace(original, '#include "stock_stub.h"')
        adapter += '\nint main(int n,char **v){if(n!=2)return 2;const char*t=getenv("TEST_THREADS");if(t)openmp_threads_num=atoi(t);if(getenv("TEST_STRESS"))stress_flag=1;Input_std(v[1]);return 0;}\n'
        (cls.build/'adapter.c').write_text(adapter)
        (cls.build/'protocol.c').write_text(PROTOCOL)
        common = [compiler, '-std=gnu99', '-I', str(SOURCE), '-I', str(UPSTREAM)]
        subprocess.run(common+[str(cls.build/'adapter.c'), str(UPSTREAM/'Inputtools.c'),
                               '-lm', '-o', str(cls.build/'adapter')], check=True,
                       capture_output=True, text=True)
        subprocess.run(common+[str(cls.build/'protocol.c'), '-o', str(cls.build/'protocol')],
                       check=True, capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='openmx run with spaces $ ')
        self.addCleanup(self.directory.cleanup)
        self.work = Path(self.directory.name)
        self.helper = self.work/'helper with spaces.py'
        self.helper.write_text(HELPER)
        self.env = {k:v for k,v in os.environ.items() if not k.startswith('OPENMX_CPU_')}
        self.env.update(OPENMX_CPU_SYMMETRY_HELPER=str(self.helper),
                        OPENMX_CPU_SYMMETRY_PYTHON=sys.executable)

    def execute(self, text, extra=None):
        input_path = self.work/'input with spaces.dat'
        input_path.write_text(text)
        env = dict(self.env)
        env.update(extra or {})
        return subprocess.run([str(self.build/'adapter'),str(input_path)], cwd=self.work,
                              env=env, capture_output=True, text=True)

    def test_environment_protocol_validation_and_profile_cleanup(self):
        result = subprocess.run([str(self.build/'protocol')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('PROTOCOL_PASS', result.stdout)

    def test_on_spawns_once_with_space_paths_and_applies_complete_profile(self):
        result = self.execute('scf.Symmetry on\n', {'OPENMX_CPU_GSOC':'1'})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('NSOC_SYM=1 GSOC=unset MODE=auto TEST=/a path=value', result.stdout)
        directories = list(self.work.glob('.openmx-symmetry-*'))
        self.assertEqual(len(directories), 1)
        self.assertEqual(directories[0].stat().st_mode & 0o777, 0o700)
        self.assertEqual((directories[0]/'bootstrap.log').read_text().count('BOOTSTRAP_CHILD_ONCE'), 1)

    def test_soc_stage_uses_spinor_profile(self):
        result = self.execute('scf.Symmetry on\nscf.SpinOrbit.Coupling on\nHS.fileout on\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('NSOC_SYM=unset GSOC=1', result.stdout)
        self.assertIn('PROJECT=1 EXPORT=1', result.stdout)

    def test_prepared_symmetry_report_reaches_native_stdout_once(self):
        result = self.execute('scf.Symmetry on\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count('CPU_SYMMETRY_SUMMARY stage=nsoc'), 1)
        self.assertEqual(result.stdout.count('CPU_SYMMETRY_OPERATION index=1/2'), 1)
        self.assertNotIn('BOOTSTRAP_CHILD_ONCE', result.stdout)

    def test_native_threads_reach_serial_preparation_helper(self):
        result = self.execute('scf.Symmetry on\n', {'TEST_THREADS':'4'})
        self.assertEqual(result.returncode, 0, result.stderr)
        log=next(self.work.glob('.openmx-symmetry-*/bootstrap.log')).read_text()
        self.assertIn('PARENT_THREADS=4',log)

    def test_verified_legacy_soc_profile_is_consumed(self):
        result=self.execute('scf.Symmetry on\nscf.SpinOrbit.Coupling on\nHS.fileout on\ntest.route legacy\n')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('GSOC=0',result.stdout)

    def test_stress_diagnostics_rejected_before_automatic_cache_selection(self):
        result=self.execute('scf.Symmetry on\n',{'TEST_STRESS':'1'})
        self.assertEqual(result.returncode,95)
        self.assertIn('stress',result.stderr.lower())

    def test_stock_logical_aliases_select_soc_stage_and_hs_export(self):
        for value in ('on','true','.true.','yes','ok','TRUE'):
            with self.subTest(value=value):
                result = self.execute('scf.Symmetry on\nscf.SpinOrbit.Coupling '+value+
                                      '\nHS.fileout '+value+'\n')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('NSOC_SYM=unset GSOC=1', result.stdout)
                self.assertIn('PROJECT=1 EXPORT=1', result.stdout)

    def test_stock_false_logical_aliases_select_nsoc_stage(self):
        for value in ('off','false','.false.','no','ng','FALSE'):
            with self.subTest(value=value):
                result = self.execute('scf.Symmetry on\nscf.SpinOrbit.Coupling '+value+'\n')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('NSOC_SYM=1 GSOC=unset', result.stdout)

    def test_invalid_native_logical_values_stop_before_preparation(self):
        for key in ('scf.SpinOrbit.Coupling', 'HS.fileout'):
            with self.subTest(key=key):
                result = self.execute('scf.Symmetry on\n'+key+' invalid\n')
                self.assertEqual(result.returncode, 95)
                self.assertNotIn('STOCK_PARSER', result.stdout)
                self.assertNotIn('CPU_INPUT_SYMMETRY_PREPARE', result.stdout)

    def test_soc_without_hs_output_keeps_projected_solver_and_disables_export(self):
        for hs_line in ('', 'HS.fileout off\n'):
            with self.subTest(hs_line=hs_line):
                result = self.execute('scf.Symmetry on\nscf.SpinOrbit.Coupling on\n'+hs_line)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('NSOC_SYM=unset GSOC=1', result.stdout)
                self.assertIn('PROJECT=1 EXPORT=0', result.stdout)

    def test_explicit_off_disables_inherited_profile_without_helper(self):
        result = self.execute('scf.Symmetry off\n', {'OPENMX_CPU_NSOC_SYM':'1',
                              'OPENMX_CPU_SYMMETRY_HELPER':'/does/not/exist'})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('NSOC_SYM=unset', result.stdout)
        self.assertEqual(list(self.work.glob('.openmx-symmetry-*')), [])

    def test_omitted_switch_preserves_legacy_developer_profile(self):
        result = self.execute('', {'OPENMX_CPU_NSOC_SYM':'1'})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('NSOC_SYM=1', result.stdout)
        self.assertEqual(list(self.work.glob('.openmx-symmetry-*')), [])

    def test_hidden_user_prepared_profile_does_not_require_bootstrap(self):
        basis = self.work/'prepared basis'
        basis.write_text('prepared fixture')
        extra = {'OPENMX_CPU_INPUT_SYMMETRY_MODE':'user',
                 'OPENMX_CPU_INPUT_SYMMETRY_GROUP':'C2',
                 'OPENMX_CPU_INPUT_SYMMETRY_TOL':'1e-5',
                 'OPENMX_CPU_INPUT_SYMMETRY_EXPORT':'raw',
                 'OPENMX_CPU_NSOC_GENERIC_COV_TOL':'1e-10',
                 'OPENMX_CPU_NSOC_SYM':'1', 'OPENMX_CPU_NSOC_BASIS':str(basis),
                 'OPENMX_CPU_SYMMETRY_HELPER':'/does/not/exist'}
        result = self.execute('scf.Symmetry.Mode User\nscf.Symmetry.PointGroup C2\n',extra)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('NSOC_SYM=1 GSOC=unset MODE=user', result.stdout)
        self.assertEqual(list(self.work.glob('.openmx-symmetry-*')), [])

    def test_on_honors_hidden_explicit_user_group(self):
        result = self.execute('scf.Symmetry on\nscf.Symmetry.Mode User\nscf.Symmetry.PointGroup C2\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('mode=User selected_group=C2', result.stdout)

    def test_invalid_toggle_and_explicit_mode_conflict_are_rejected(self):
        for text in ('scf.Symmetry garbage\n',
                     'scf.Symmetry off\nscf.Symmetry.Mode Auto\n',
                     'scf.Symmetry on\nscf.Symmetry.Mode Off\n'):
            with self.subTest(text=text):
                result = self.execute(text)
                self.assertEqual(result.returncode, 95)
                self.assertNotIn('STOCK_PARSER', result.stdout)

    def test_child_failure_is_rejected_before_stock_parser(self):
        result = self.execute('scf.Symmetry on\ntest.helper fail\n')
        self.assertEqual(result.returncode, 95)
        self.assertIn('preparation failed', result.stderr)
        self.assertNotIn('STOCK_PARSER', result.stdout)

    def test_oversized_environment_is_rejected_before_stock_parser(self):
        result = self.execute('scf.Symmetry on\ntest.helper oversize\n')
        self.assertEqual(result.returncode, 95)
        self.assertIn('invalid preparation environment', result.stderr)
        self.assertNotIn('STOCK_PARSER', result.stdout)


if __name__ == '__main__':
    unittest.main()
