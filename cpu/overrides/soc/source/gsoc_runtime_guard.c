/* Collective validation before any GSOC-dependent allocation or collective.
 * Digests detect rank disagreement; they are not authentication primitives. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include "gsoc_runtime_guard.h"
#ifndef GSOC_RUNTIME_GUARD_STANDALONE
#include "openmx_common.h"
extern int flag_export_DM,flag_energy_range_DM;
#endif
static void digest_bytes(uint64_t h[2],const void *data,size_t size){
 const unsigned char *p=(const unsigned char*)data;
 for(size_t i=0;i<size;i++){h[0]=(h[0]^p[i])*UINT64_C(1099511628211);h[1]=(h[1]+p[i]+UINT64_C(0x9e3779b97f4a7c15))*UINT64_C(0xbf58476d1ce4e5b9);h[1]^=h[1]>>29;}
}
static void digest_init(uint64_t h[2]){h[0]=UINT64_C(14695981039346656037);h[1]=UINT64_C(0x94d049bb133111eb);}
static int error_text(char *out,size_t size,const char *text){if(size)snprintf(out,size,"%s",text);return -1;}
static int int_option(const char *key,int fallback,int low,int high,int *out,char *error,size_t size){
 const char *s=getenv(key);if(!s){*out=fallback;return 0;}
 if(low==0&&high==1&&strcmp(s,"0")&&strcmp(s,"1"))return error_text(error,size,key);
 char *end;errno=0;long x=strtol(s,&end,10);
 if(end==s||*end||errno||x<low||x>high)return error_text(error,size,key);
 *out=(int)x;return 0;
}
static int real_option(const char *key,double fallback,double low,double high,double *out,char *error,size_t size){
 const char *s=getenv(key);if(!s){*out=fallback;return 0;}
 char *end;errno=0;double x=strtod(s,&end);
 if(end==s||*end||errno||!isfinite(x)||x<low||x>high)return error_text(error,size,key);
 *out=x;return 0;
}
/* Hash the explicit communication/lifetime options only, without exposing
 * values or requiring identical environment ordering. Basis paths may differ
 * when their complete contents agree. */
static void digest_environment(uint64_t h[2]){
 const char *keys[]={"OPENMX_CPU_GSOC","OPENMX_CPU_GSOC_COMPACT",
  "OPENMX_CPU_GSOC_TILE","OPENMX_CPU_GSOC_PROBES","OPENMX_CPU_GSOC_PROJECT",
  "OPENMX_CPU_GSOC_ALLOW_FRESH","OPENMX_CPU_GSOC_BLOCK_BUCKETS",
  "OPENMX_CPU_GSOC_COV_TOL","OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG",
  "OPENMX_CPU_C3","OPENMX_CPU_D3","OPENMX_CPU_SOC_COMPACT",
  "OPENMX_CPU_SOC_TR","OPENMX_CPU_SOC_REAL","OPENMX_CPU_SOC_RELEASE_DENSE",
  "OPENMX_CPU_DM_BLAS3","OPENMX_CPU_DM_SORT","OPENMX_CPU_DM_TRIM",
  "OPENMX_CPU_DM_CYCLIC","OPENMX_CPU_DM_BALANCE_PROFILE"};
 digest_init(h);
 for(size_t k=0;k<sizeof(keys)/sizeof(*keys);k++){
  const char *value=getenv(keys[k]);unsigned char present=value?1:0;
  digest_bytes(h,&present,1);if(value)digest_bytes(h,value,strlen(value)+1);
 }
}
static int read_basis(GsocRuntimeConfig *c,const char *path,char *error,size_t size){
 if(!path||!*path)return error_text(error,size,"missing OPENMX_CPU_GSOC_BASIS");
 FILE *f=fopen(path,"rb");if(!f)return error_text(error,size,"cannot read GSOC basis");
 int status=-1;
 if(fread(c->header,4,6,f)!=6||fread(&c->basis_nnz,8,1,f)!=1)goto bad_header;
 int32_t *h=c->header;
 if(h[0]!=0x47534f43||(h[1]!=1&&h[1]!=2)||h[2]<1||h[3]<1||h[3]>INT_MAX/4||h[4]<1||h[4]>48||h[5]<1||h[5]>192||c->basis_nnz<1||c->basis_nnz>INT_MAX)goto bad_header;
 if(c->input.atoms>0&&h[2]!=c->input.atoms)goto bad_header;
 if(fread(c->sector_size,4,h[4],f)!=(size_t)h[4]||fread(c->sector_type,4,h[4],f)!=(size_t)h[4]||fread(c->sector_partner,4,h[4],f)!=(size_t)h[4])goto bad_header;
 int64_t states=0;
 for(int k=0;k<h[4];k++){
  int n=c->sector_size[k],p=c->sector_partner[k],t=c->sector_type[k];
  if(n<1||(t!=1&&t!=2)||p< -1||p>=k)goto bad_header;
  if(p>=0&&(c->sector_partner[p]>=0||n!=c->sector_size[p]||t!=c->sector_type[p]))goto bad_header;
  if(p>=0)for(int j=0;j<k;j++)if(c->sector_partner[j]==p)goto bad_header;
  states+=n;
 }
 if(states!=2*(int64_t)h[3])goto bad_header;
 if(fseek(f,0,SEEK_SET))goto read_failed;
 digest_init(c->basis_digest);unsigned char buffer[65536];size_t n;
 while((n=fread(buffer,1,sizeof(buffer),f))>0){
  if(c->basis_bytes>UINT64_MAX-n)goto read_failed;
  c->basis_bytes+=n;digest_bytes(c->basis_digest,buffer,n);
 }
 if(ferror(f))goto read_failed;
 /* Require the complete declared payload, not just a plausible sector header. */
 uint64_t expected=32+12*(uint64_t)h[4]+72+8*(uint64_t)h[2]
  +(uint64_t)h[5]*(96+4*(uint64_t)h[2])+8*(2*(uint64_t)h[3]+1)+20*(uint64_t)c->basis_nnz;
 if(c->basis_bytes!=expected)goto bad_header;
 status=0;goto done;
 bad_header:error_text(error,size,"invalid or incomplete GSOC basis/sector header");goto done;
 read_failed:error_text(error,size,"failed to hash complete GSOC basis");
 done:fclose(f);return status;
}
int gsoc_runtime_capture(GsocRuntimeConfig *c,const GsocRuntimeInput *input,char *error,size_t size){
 memset(c,0,sizeof(*c));if(size)*error=0;
 if(int_option("OPENMX_CPU_GSOC",0,0,1,&c->enabled,error,size)
  ||int_option("OPENMX_CPU_GSOC_COMPACT",0,0,1,&c->compact,error,size))return -1;
 if(c->compact&&!c->enabled)return error_text(error,size,"GSOC_COMPACT requires GSOC=1; disabled raw path needs full EVec");
 /* Off must not require a basis or impose GSOC restrictions on ordinary runs. */
 if(!c->enabled)return 0;
 c->input=*input;
 int legacy=0;
 const char *legacy_keys[]={"OPENMX_CPU_C3","OPENMX_CPU_D3","OPENMX_CPU_SOC_COMPACT","OPENMX_CPU_SOC_TR","OPENMX_CPU_SOC_REAL"};
 for(size_t k=0;k<sizeof(legacy_keys)/sizeof(*legacy_keys);k++){
  if(int_option(legacy_keys[k],0,0,k==0?6:1,&legacy,error,size))return -1;
  if(legacy)return error_text(error,size,"GSOC requires legacy C3/D3/compact/TR/real selectors off");
 }
 if(int_option("OPENMX_CPU_GSOC_TILE",32,1,128,&c->tile,error,size)
  ||int_option("OPENMX_CPU_GSOC_PROBES",4,1,16,&c->probes,error,size)
  ||int_option("OPENMX_CPU_GSOC_PROJECT",0,0,1,&c->project,error,size)
  ||int_option("OPENMX_CPU_GSOC_ALLOW_FRESH",0,0,1,&c->allow_fresh,error,size)
  ||int_option("OPENMX_CPU_GSOC_BLOCK_BUCKETS",1,0,1,&c->block_buckets,error,size)
  ||real_option("OPENMX_CPU_GSOC_COV_TOL",1e-10,1e-14,1e-3,&c->covariance_tol,error,size)
  ||real_option("OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG",5e-5,1e-12,0.02,&c->geometry_tol,error,size))return -1;
 if(!c->project&&c->covariance_tol>1e-10)return error_text(error,size,"relaxed covariance requires explicit GSOC_PROJECT=1");
 if(input->solver!=2||input->spin!=3||input->soc!=1||input->md||input->contracted||input->scf_steps!=1)return error_text(error,size,"GSOC requires fixed one-step Gamma cluster SOC and uncontracted basis");
 if(input->hubbard||input->constraint||input->zeeman_spin||input->zeeman_orbital)return error_text(error,size,"GSOC excludes Hubbard, constraint and Zeeman spin extensions");
 if(input->c2n&&(input->restart!=1||input->restart_spin!=0))return error_text(error,size,"GSOC requires nonmagnetic scalar c2n checkpoint");
 if(!input->c2n&&!c->allow_fresh)return error_text(error,size,"GSOC requires c2n or explicit small fresh smoke flag");
 if(input->dense_postprocessing)return error_text(error,size,"GSOC does not support requested native state-based postprocessing");
 if(input->cwf||(input->export_dm>0&&input->energy_window_dm>0)||input->dos||input->dos_gauss||input->optical)return error_text(error,size,"GSOC dense placeholders forbid CWF, energy-window DM, DOS re-solves and optical postprocessing");
 if(input->partial_charge||input->atomic_analysis)return error_text(error,size,"GSOC excludes partial-charge and atomic-orbital analysis postprocessing");
 if(c->compact){
  const char *dm=getenv("OPENMX_CPU_DM_BLAS3");
  if(!dm||strcmp(dm,"1"))return error_text(error,size,"compact GSOC requires DM_BLAS3 exactly1; audit mode2 needs full EVec");
  if(input->level_fileout>=2)return error_text(error,size,"compact GSOC does not provide full eigenvector printout");
  const char *balance_keys[]={"OPENMX_CPU_DM_SORT","OPENMX_CPU_DM_TRIM","OPENMX_CPU_DM_CYCLIC"};
  for(size_t k=0;k<sizeof(balance_keys)/sizeof(*balance_keys);k++){
   int value=0;if(int_option(balance_keys[k],0,0,1,&value,error,size))return -1;
   if(value)return error_text(error,size,"initial compact GSOC requires sort/trim/cyclic off");
  }
 }
 if(c->compact&&input->mo)return error_text(error,size,"compact GSOC does not provide full MO coefficient export");
 digest_environment(c->environment_digest);
 if(read_basis(c,getenv("OPENMX_CPU_GSOC_BASIS"),error,size))return -1;
 if(!input->c2n&&2*c->header[3]>512)return error_text(error,size,"fresh GSOC is restricted to <=512 spinors");
 return 0;
}
void gsoc_runtime_agree(MPI_Comm comm,const GsocRuntimeConfig *c,int local_status,const char *error){
 int rank,failed=local_status?1:0,any_failed=0;MPI_Comm_rank(comm,&rank);
 MPI_Allreduce(&failed,&any_failed,1,MPI_INT,MPI_MAX,comm);
 if(any_failed){if(failed){fprintf(stderr,"CPU_GSOC_RUNTIME_ERROR rank=%d %s\n",rank,error?error:"invalid local configuration");fflush(stderr);}MPI_Abort(comm,95);abort();}
 uint64_t local[2],reference[2];digest_init(local);digest_bytes(local,c,sizeof(*c));
 if(!rank){reference[0]=local[0];reference[1]=local[1];}
 MPI_Bcast(reference,2,MPI_UINT64_T,0,comm);
 int mismatch=local[0]!=reference[0]||local[1]!=reference[1],any_mismatch=0;
 MPI_Allreduce(&mismatch,&any_mismatch,1,MPI_INT,MPI_MAX,comm);
 if(any_mismatch){if(mismatch){fprintf(stderr,"CPU_GSOC_RUNTIME_ERROR rank=%d MPI configuration/basis disagreement\n",rank);fflush(stderr);}MPI_Abort(comm,95);abort();}
 if(c->enabled&&!rank){printf("CPU_GSOC_RUNTIME_GUARD ranks_agree=1 compact=%d basis_bytes=%llu sectors=%d tile=%d probes=%d\n",c->compact,(unsigned long long)c->basis_bytes,c->header[4],c->tile,c->probes);fflush(stdout);}
}
#ifndef GSOC_RUNTIME_GUARD_STANDALONE
void CpuGSOC_ValidateRuntime(void){
 GsocRuntimeInput in;GsocRuntimeConfig c;char error[256];memset(&in,0,sizeof(in));
 in.solver=Solver;in.spin=SpinP_switch;in.soc=SO_switch;in.md=MD_switch;
 in.contracted=Cnt_switch;in.scf_steps=DFTSCF_loop;in.c2n=Use_of_Collinear_Restart;
 in.restart=Scf_RestartFromFile;in.restart_spin=SpinP_switch_RestartFiles;
 in.hubbard=Hub_U_switch;in.constraint=Constraint_NCS_switch;
 in.zeeman_spin=Zeeman_NCS_switch;in.zeeman_orbital=Zeeman_NCO_switch;
 in.cwf=CWF_Calc;in.export_dm=flag_export_DM;in.energy_window_dm=flag_energy_range_DM;
 in.level_fileout=level_fileout;in.partial_charge=cal_partial_charge;
 in.atomic_analysis=LNAO_calc_flag||LNBO_calc_flag||NBO_switch;
 in.dense_postprocessing=Band_disp_switch||Wannier_Func_Calc||Wannier90_fileout||CWF_fileout_flag
  ||pop_anal_aow_flag||unfold_electronic_band||Polarization_ABC_flag||xanes_gs_fileout||xanes_calc;
 in.dos=Dos_fileout;in.dos_gauss=DosGauss_fileout;in.optical=CDDF_on;in.mo=MO_fileout;in.atoms=atomnum;
 int status=gsoc_runtime_capture(&c,&in,error,sizeof(error));
 gsoc_runtime_agree(mpi_comm_level1,&c,status,error);
}
#endif
