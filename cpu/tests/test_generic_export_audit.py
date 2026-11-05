"""Compile and execute the independent solved-sector export matvec guard."""
from pathlib import Path
import shutil,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
class ExportAudit(unittest.TestCase):
 def test_sector_TR_pair_and_real_projection_match_dense_oracle(self):
  code=r'''
#include "gsoc_export_audit.h"
static GsocAuditComplex H[16],S[16];
static int basis(void*c,const GsocAuditComplex*x,GsocAuditComplex*y){memcpy(y,x,4*sizeof(*x));return 0;}
static int apply(void*c,const GsocAuditComplex*x,GsocAuditComplex*h,GsocAuditComplex*s){
 for(int i=0;i<4;i++){h[i]=s[i]=gsoc_audit_c(0,0);for(int j=0;j<4;j++){
  GsocAuditComplex a=gsoc_audit_mul(H[4*i+j],x[j]),b=gsoc_audit_mul(S[4*i+j],x[j]);h[i].r+=a.r;h[i].i+=a.i;s[i].r+=b.r;s[i].i+=b.i;}}return 0;}
int main(){int ofs[3]={0,2,4},dtype[2]={1,1},partner[2]={-1,0};
 GsocAuditOperator o={4,2,ofs,dtype,0,apply,basis,basis};
 GsocAuditComplex x[4]={{.2,.3},{-.7,.8},{.4,-.2},{.5,.6}},h[4],s[4];
 for(int i=0;i<4;i++)for(int j=0;j<4;j++){H[4*i+j]=gsoc_audit_c(1+i+j,(i-j)*.17);S[4*i+j]=gsoc_audit_c(i==j?2:.1,0);}
 if(gsoc_export_projected_apply(&o,partner,x,h,s))return 1;
 for(int i=0;i<4;i++){GsocAuditComplex want={0,0},sw={0,0};for(int j=0;j<2;j++){
  int ii=i%2;double a=.5*(H[4*ii+j].r+H[4*(ii+2)+j+2].r),b=S[4*ii+j].r;
  want.r+=a*x[2*(i/2)+j].r;want.i+=a*x[2*(i/2)+j].i;sw.r+=b*x[2*(i/2)+j].r;sw.i+=b*x[2*(i/2)+j].i;}
  if(hypot(h[i].r-want.r,h[i].i-want.i)>1e-12||hypot(s[i].r-sw.r,s[i].i-sw.i)>1e-12)return 2;}
 dtype[0]=dtype[1]=2;partner[1]=-1;
 if(gsoc_export_projected_apply(&o,partner,x,h,s))return 3;
 for(int i=0;i<4;i++){GsocAuditComplex want={0,0};for(int j=2*(i/2);j<2*(i/2)+2;j++){
  GsocAuditComplex a=gsoc_audit_mul(H[4*i+j],x[j]);want.r+=a.r;want.i+=a.i;}
  if(hypot(h[i].r-want.r,h[i].i-want.i)>1e-12)return 4;}
 return 0;
}'''
  with tempfile.TemporaryDirectory()as d:
   source=Path(d)/'test.c';source.write_text(code);exe=Path(d)/'test'
   subprocess.run([shutil.which('cc')or'cc','-std=c99','-I'+str(ROOT/'overrides/soc/source'),str(source),'-lm','-o',str(exe)],check=True,capture_output=True,text=True)
   subprocess.run([str(exe)],check=True)
if __name__=='__main__':unittest.main()
