"""Numerical regression: raw symmetry defects versus a projected eigensolve."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

CPU = Path(__file__).resolve().parents[1]


class ProjectedSocPolicy(unittest.TestCase):
    def test_projected_operator_keeps_solve_checks_without_raw_spatial_gate(self):
        code = r'''
#include <assert.h>
#include "gsoc_export_audit.h"
static GsocAuditComplex H[16], S[16];
static int basis(void*c,const GsocAuditComplex*x,GsocAuditComplex*y){memcpy(y,x,4*sizeof(*x));return 0;}
static int apply(void*c,const GsocAuditComplex*x,GsocAuditComplex*h,GsocAuditComplex*s){
 for(int i=0;i<4;i++){h[i]=s[i]=gsoc_audit_c(0,0);for(int j=0;j<4;j++){
  GsocAuditComplex a=gsoc_audit_mul(H[4*i+j],x[j]),b=gsoc_audit_mul(S[4*i+j],x[j]);
  h[i].r+=a.r;h[i].i+=a.i;s[i].r+=b.r;s[i].i+=b.i;}}return 0;}
static int bad_eigenvalue=0;
static int vector(void*c,int band,GsocAuditComplex*x,double*e){
 memset(x,0,4*sizeof(*x));x[band-1].r=1;*e=(band-1)%2+1+bad_eigenvalue;return 0;}
int main(void){
 int ofs[3]={0,2,4},dtype[2]={1,1},partner[2]={-1,0},bands[4]={1,2,3,4};
 GsocAuditOperator raw={4,2,ofs,dtype,0,apply,basis,basis};
 for(int i=0;i<4;i++){H[4*i+i]=gsoc_audit_c(i%2+1,0);S[4*i+i]=gsoc_audit_c(1,0);}
 /* Real antisymmetric up/down coupling preserves physical J K but crosses sectors. */
 H[3]=H[12]=gsoc_audit_c(.1,0);H[6]=H[9]=gsoc_audit_c(-.1,0);
 S[3]=S[12]=gsoc_audit_c(.02,0);S[6]=S[9]=gsoc_audit_c(-.02,0);
 GsocAuditGateResult gate;assert(!gsoc_audit_operator_gate(&raw,4,20261106,&gate));
 assert(gate.comm_rel[0]>1e-3 && gate.comm_rel[1]>1e-3);
 assert(gate.tr_rel[0]<1e-12 && gate.tr_rel[1]<1e-12);
 assert(gsoc_audit_gate_status(&gate,0,1e-5)!=0);
 assert(gsoc_audit_gate_status(&gate,1,1e-5)==0);
 gate.real_rel[0]=.1;assert(gsoc_audit_gate_status(&gate,1,1e-5)==0);
 assert(gsoc_audit_gate_status(&gate,0,1e-5)!=0);gate.real_rel[0]=0;
 gate.tr_rel[0]=.01;assert(gsoc_audit_gate_status(&gate,1,1e-5)!=0);gate.tr_rel[0]=0;
 gate.nonfinite=1;assert(gsoc_audit_gate_status(&gate,1,1e-5)!=0);gate.nonfinite=0;
 gate.comm_rel[0]=NAN;assert(gsoc_audit_gate_status(&gate,1,1e-5)!=0);
 GsocProjectedAudit context={&raw,partner};GsocAuditOperator solved=raw;
 solved.ctx=&context;solved.apply=gsoc_export_projected_callback;
 GsocAuditEigenResult a,b;assert(!gsoc_audit_eigenpairs(&raw,vector,bands,4,&a));
 assert(a.raw_relative>1e-3 && a.Sorth>1e-3);
 assert(!gsoc_audit_eigenpairs(&solved,vector,bands,4,&b));
 assert(b.raw_relative<1e-12 && b.Sorth<1e-12 && !b.nonfinite);
 bad_eigenvalue=1;assert(!gsoc_audit_eigenpairs(&solved,vector,bands,4,&b));
 assert(b.raw_relative>1e-3);return 0;
}'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'test.c').write_text(code)
            result = subprocess.run([shutil.which('cc') or 'cc', '-std=c99',
                '-I'+str(CPU/'overrides/soc/source'), str(root/'test.c'), '-lm',
                '-o', str(root/'test')], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            subprocess.run([str(root/'test')], check=True)


if __name__ == '__main__':
    unittest.main()
