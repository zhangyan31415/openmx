"""Compile the shared periodic gauge rebase without MPI or a native SCF."""
from pathlib import Path
import subprocess,tempfile,unittest,shutil
CPU=Path(__file__).resolve().parents[1]
C=r'''
#include <assert.h>
#include <string.h>
#include "cpu_periodic_image_gauge.h"
int main(void){
 double cell[9]={4,0,0,1,5,0,.5,.2,6};
 double native_data[4][4]={{0},{0,.4,.6,.8},{0,1.3,1.5,1.7},{0,2.1,2.3,2.5}};
 double *actual[4]={native_data[0],native_data[1],native_data[2],native_data[3]};
 int m[3][3]={{2,-1,3},{-4,2,-1},{1,0,-2}};
 double ref[9];for(int a=0;a<3;a++)for(int i=0;i<3;i++){ref[3*a+i]=actual[a+1][i+1];for(int j=0;j<3;j++)ref[3*a+i]+=m[a][j]*cell[3*j+i];}
 int p0[3]={0,1,2},p1[3]={2,1,0};int *perm[2]={p0,p1};
 int w[2][9]={{1,0,0,0,1,0,0,0,1},{-1,0,0,0,-1,0,0,0,-1}};
 int s0[9]={0},s1[9]={1,2,3,4,5,6,7,8,9},old[9];memcpy(old,s1,sizeof(old));int *image[2]={s0,s1};char error[256];int moved=-1;
 assert(!cpu_periodic_image_rebase(3,2,cell,ref,actual,1.,perm,w,image,2e-7,&moved,error,sizeof(error)));assert(moved==3);
 for(int a=0;a<3;a++)for(int i=0;i<3;i++){
  assert(ref[3*a+i]==actual[a+1][i+1]);assert(s0[3*a+i]==0);
  assert(s1[3*a+i]==old[3*a+i]+m[p1[a]][i]+m[a][i]);
 }
 memcpy(old,s1,sizeof(old));assert(!cpu_periodic_image_rebase(3,2,cell,ref,actual,1.,perm,w,image,2e-7,&moved,error,sizeof(error)));assert(!moved);assert(!memcmp(old,s1,sizeof(old)));
 ref[0]+=.01;assert(cpu_periodic_image_rebase(3,2,cell,ref,actual,1.,perm,w,image,2e-7,&moved,error,sizeof(error)));assert(strstr(error,"not a lattice translation"));
 return 0;
}
'''
class Gauge(unittest.TestCase):
 def test_image_gauge_rebases_arbitrary_periodic_representatives_and_rejects_real_motion(self):
  header=CPU/'overrides/soc/source/cpu_periodic_image_gauge.h';self.assertTrue(header.exists())
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);(root/'test.c').write_text(C);cc=shutil.which('cc')or shutil.which('gcc');self.assertIsNotNone(cc)
   subprocess.run([cc,'-std=c99','-O2','-Wall','-Werror','-I',str(header.parent),str(root/'test.c'),'-lm','-o',str(root/'test')],check=True)
   subprocess.run([str(root/'test')],check=True)
if __name__=='__main__':unittest.main()
