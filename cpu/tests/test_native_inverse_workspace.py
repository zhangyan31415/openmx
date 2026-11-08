"""Regression: Pulay LU scratch must not use the physical DM slot count."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'overrides/nsoc/source/DIIS_Mixing_DM.unified.c'


class NativeInverseWorkspace(unittest.TestCase):
    @unittest.skipUnless(shutil.which('cc'), 'C compiler unavailable')
    def test_full_history_with_three_physical_dm_slots(self):
        source = SOURCE.read_text()
        start = source.index('\nvoid Inverse(int n, double **a, double **ia)\n{')
        end = source.index('\nvoid DIIS_Mixing_DM(', start)
        fixture = r'''
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <sys/mman.h>
#include <unistd.h>
/* A guard directly after every allocation makes even a one-element overrun
 * fail without relying on optional ASan libraries or allocator behavior. */
static void *guarded_malloc(size_t n){
 size_t pg=(size_t)sysconf(_SC_PAGESIZE);assert(n<=pg);
 char *p=mmap(NULL,pg*2,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
 assert(p!=MAP_FAILED);assert(!mprotect(p+pg,pg,PROT_NONE));return p+pg-n;
}
static void guarded_free(void *p){(void)p;}
#define malloc guarded_malloc
#define free guarded_free
int List_YOUSO[100];
''' + source[start:end] + r'''
int main(void){
 List_YOUSO[16]=3;
 for(int n=1;n<=30;n++){
  double data[32][32]={{0}}, original[32][32],inv[32][32],*a[32],*ia[32];
  for(int i=0;i<=n;i++){
   a[i]=data[i];ia[i]=inv[i];
   for(int j=0;j<=n;j++){
    data[i][j]=(i==n&&j==n)?0.:(i==n||j==n)?-1.:(i==j)?i+2.:0.;
    original[i][j]=data[i][j];
   }
  }
  Inverse(n,a,ia);
  for(int i=0;i<=n;i++)for(int j=0;j<=n;j++){
   double v=0.;for(int k=0;k<=n;k++)v+=original[i][k]*inv[k][j];
   assert(fabs(v-(i==j))<1e-10);assert(data[i][j]==original[i][j]);
  }
 }
 return 0;
}
'''
        with tempfile.TemporaryDirectory(prefix='openmx-inverse-') as directory:
            root = Path(directory)
            (root/'inverse.c').write_text(fixture)
            subprocess.run(['cc', '-O1', str(root/'inverse.c'), '-o', str(root/'inverse')], check=True, capture_output=True)
            subprocess.run([str(root/'inverse')], check=True, capture_output=True)


if __name__ == '__main__':
    unittest.main()
