"""Collective preflight captures and hashes full generic export maps."""
from pathlib import Path
import os,shutil,struct,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
class ExportGuard(unittest.TestCase):
 def test_capture_rejects_missing_projection_and_detects_map_changes(self):
  with tempfile.TemporaryDirectory()as d:
   d=Path(d)
   (d/'mpi.h').write_text('typedef int MPI_Comm;\n#define MPI_INT 1\n#define MPI_MAX 2\n#define MPI_UINT64_T 3\nint MPI_Comm_rank();int MPI_Allreduce();int MPI_Bcast();int MPI_Abort();\n')
   (d/'test.c').write_text('''#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gsoc_runtime_guard.h"
int main(){GsocRuntimeConfig c;GsocRuntimeInput i;char error[256];memset(&i,0,sizeof(i));i.solver=2;i.spin=3;i.soc=1;i.scf_steps=1;i.c2n=1;i.restart=1;i.atoms=1;int s=gsoc_runtime_capture(&c,&i,error,sizeof(error));if(s){puts(error);return 1;}printf("%llu %llu\\n",(unsigned long long)c.map_digest[0],(unsigned long long)c.map_digest[1]);return 0;}
''')
   exe=d/'test';subprocess.run([shutil.which('cc')or'cc','-std=c99','-DGSOC_RUNTIME_GUARD_STANDALONE','-ffunction-sections','-I'+str(d),'-I'+str(ROOT/'overrides/soc/source'),str(d/'test.c'),str(ROOT/'overrides/soc/source/gsoc_runtime_guard.c'),'-Wl,--gc-sections','-lm','-o',str(exe)],check=True,capture_output=True,text=True)
   basis=d/'basis';payload=struct.pack('<6iq',0x47534f43,1,1,1,2,1,2)+struct.pack('<6i',1,1,2,2,-1,-1)
   # Header/sector arrays followed by complete declared opaque payload.
   basis.write_bytes(payload+bytes(300-len(payload)))
   meta=d/'map';header=struct.pack('<8i',0x52534433,3,1,1,1,1,0,0)
   payload=header+bytes(72+28)+struct.pack('<i',1)+bytes(36+16+8+64);self.assertEqual(len(payload),260);meta.write_bytes(payload)
   env={k:v for k,v in os.environ.items()if not k.startswith('OPENMX_CPU_')};env.update(OPENMX_CPU_GSOC='1',OPENMX_CPU_GSOC_BASIS=str(basis),OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE='1',OPENMX_CPU_SOC_NATIVE_MAP=str(meta))
   r=subprocess.run([str(exe)],env=env,capture_output=True,text=True);self.assertNotEqual(r.returncode,0);self.assertIn('PROJECT=1',r.stdout)
   env['OPENMX_CPU_GSOC_PROJECT']='1';r=subprocess.run([str(exe)],env=env,capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stdout)
   changed=bytearray(payload);changed[-1]=1;meta.write_bytes(changed)
   q=subprocess.run([str(exe)],env=env,capture_output=True,text=True);self.assertEqual(q.returncode,0,q.stdout);self.assertNotEqual(q.stdout,r.stdout)
   meta.write_bytes(changed[:-1]);q=subprocess.run([str(exe)],env=env,capture_output=True,text=True);self.assertNotEqual(q.returncode,0);self.assertIn('incomplete',q.stdout)
if __name__=='__main__':unittest.main()
