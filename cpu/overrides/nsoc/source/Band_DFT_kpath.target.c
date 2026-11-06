/* Targeted band dispersion: retain the stock full-spectrum implementation. */
#define Band_DFT_kpath CpuOriginal_Band_DFT_kpath
#include "../../../vendor/openmx4.0/source/Band_DFT_kpath.c"
#undef Band_DFT_kpath
#include <stdint.h>
#include <limits.h>
#include "Inputtools.h"
#include "cpu_band_options.h"
#include "cpu_band_bridge.h"
#include "cpu_soc_native.h"

static void bt_fail(const char *message){
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 fprintf(stderr,"CPU_TARGET_BANDS_ERROR rank=%d %s\n",rank,message);fflush(stderr);
 MPI_Abort(mpi_comm_level1,96);abort();
}
static void *bt_mem(size_t n,size_t bytes){
 if(bytes&&n>SIZE_MAX/bytes)bt_fail("allocation overflow");
 void *p=calloc(n?n:1,bytes);if(!p)bt_fail("allocation failed");return p;
}
static void bt_write(const void *data,size_t size,size_t count,FILE *fp){if(fwrite(data,size,count,fp)!=count)bt_fail("snapshot write failed");}
static void bt_read(void *data,size_t size,size_t count,FILE *fp){if(fread(data,size,count,fp)!=count)bt_fail("short worker result");}
static void bt_parts(const char *directory,int rank,int n,int spin,double *****h,double *****im,double ****s){
 char path[PATH_MAX];snprintf(path,sizeof(path),"%s/part-%06d.bin",directory,rank);
 FILE *fp=fopen(path,"wb");if(!fp)bt_fail("cannot write native H/S snapshot");
 int *offset=bt_mem(atomnum+1,sizeof(int)),start=0,blocks=0;
 for(int a=1;a<=atomnum;a++){offset[a]=start;int na=Spe_Total_CNO[WhatSpecies[a]];if(na>INT_MAX-start)bt_fail("AO count overflow");start+=na;}
 for(int a=1;a<=Matomnum;a++){int ga=M2G[a];if(FNAN[ga]>INT_MAX-blocks-1)bt_fail("neighbor count overflow");blocks+=FNAN[ga]+1;}
 int fields=spin==3?8:(spin==1?3:2);int32_t header[5]={n,spin,blocks,fields,1};
 bt_write("OMBAND01",1,8,fp);bt_write(header,4,5,fp);
 int imaginary=SO_switch||Hub_U_switch||Constraint_NCS_switch||Zeeman_NCS_switch||Zeeman_NCO_switch;
 for(int a=1;a<=Matomnum;a++){
  int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];
  for(int b=0;b<=FNAN[ga];b++){
   int gb=natn[ga][b],nb=Spe_Total_CNO[WhatSpecies[gb]],image=ncn[ga][b];
   int32_t block[7]={offset[ga],offset[gb],na,nb,atv_ijk[image][1],atv_ijk[image][2],atv_ijk[image][3]};bt_write(block,4,7,fp);
   for(int field=0;field<fields;field++)for(int i=0;i<na;i++)for(int j=0;j<nb;j++){
    double value;
    if(field==fields-1)value=s[a][b][i][j];
    else if(spin!=3||field<4)value=h[field][a][b][i][j];
    else value=imaginary?im[field-4][a][b][i][j]:0.;
    if(!isfinite(value))bt_fail("nonfinite native H/S value");bt_write(&value,8,1,fp);
   }
  }
 }
 if(fclose(fp))bt_fail("cannot close native snapshot");free(offset);
}
static void bt_write_band(int nkpath,int *n_perk,double ***kpath,char ***kname,int channels,int count,int points,const double *values){
 char path[YOUSO10];snprintf(path,sizeof(path),"%s%s.Band",filepath,filename);FILE *fp=fopen(path,"w");
 if(!fp)bt_fail("cannot write .Band output");
 fprintf(fp," %d %d %18.15f\n",count,channels-1,ChemP);
 for(int i=1;i<=3;i++)for(int j=1;j<=3;j++)fprintf(fp,"%18.15f ",rtv[i][j]);fprintf(fp,"\n%d\n",nkpath);
 for(int p=1;p<=nkpath;p++)fprintf(fp,"%d %18.15f %18.15f %18.15f %18.15f %18.15f %18.15f %s %s\n",n_perk[p],kpath[p][1][1],kpath[p][1][2],kpath[p][1][3],kpath[p][2][1],kpath[p][2][2],kpath[p][2][3],kname[p][1],kname[p][2]);
 int point=0;
 for(int p=1;p<=nkpath;p++)for(int j=0;j<n_perk[p];j++,point++)for(int channel=0;channel<channels;channel++){
  fprintf(fp,"%d",count);
  for(int d=1;d<=3;d++)fprintf(fp," %18.15f",kpath[p][1][d]+(kpath[p][2][d]-kpath[p][1][d])*(double)j/(n_perk[p]-1));
  fprintf(fp,"\n");for(int b=0;b<count;b++)fprintf(fp,"%18.15f ",values[((size_t)point*channels+channel)*count+b]);fprintf(fp,"\n");
 }
 if(point!=points||fclose(fp))bt_fail("band output count/close failed");
}
static void bt_run_bands(int nkpath,int *n_perk,double ***kpath,char ***kname,int spin,double *****h,double *****im,double ****s){
 if(!CpuBandOptions.count){
  double start=MPI_Wtime();CpuOriginal_Band_DFT_kpath(nkpath,n_perk,kpath,kname,spin,h,im,s);
  double local=MPI_Wtime()-start,maximum;int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
  MPI_Reduce(&local,&maximum,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_level1);
  if(!rank){printf("CPU_FULL_BANDS_DONE total_seconds=%.6f full_spectrum=1\n",maximum);fflush(stdout);}return;
 }
 int rank,np;MPI_Comm_rank(mpi_comm_level1,&rank);MPI_Comm_size(mpi_comm_level1,&np);
 if(Cnt_switch||LNO_flag||!(spin==0||spin==1||spin==3)||!h||!s)bt_fail("targeted native bands require an uncontracted ordinary AO basis (no LNO)");
 double started=MPI_Wtime();int n=0,points=0,channels=spin==1?2:1,count=CpuBandOptions.count;
 for(int a=1;a<=atomnum;a++){int na=Spe_Total_CNO[WhatSpecies[a]];if(na>INT_MAX-n)bt_fail("AO count overflow");n+=na;}
 if(n>INT_MAX/2)bt_fail("spinor dimension overflow");
 int dimension=spin==3?2*n:n;if(count>=dimension-1)bt_fail("Band.NumBands must be smaller than matrix dimension minus one");
 if(spin==3&&(SO_switch||Hub_U_switch||Constraint_NCS_switch||Zeeman_NCS_switch||Zeeman_NCO_switch)&&!im)bt_fail("missing native SOC imaginary terms");
 for(int p=1;p<=nkpath;p++){if(n_perk[p]<2||n_perk[p]>INT_MAX-points)bt_fail("each Band.kpath segment needs at least two points");points+=n_perk[p];}
 if(!points||points>INT_MAX/channels/count)bt_fail("band result MPI count overflow");
 int workers=points<np?points:np;double target=ChemP;
 if(strcasecmp(CpuBandOptions.reference,"Fermi"))target=strtod(CpuBandOptions.reference,NULL)/27.211386245988;
 if(!isfinite(target))bt_fail("nonfinite reference energy");
 char directory[PATH_MAX]={0};
 if(!rank){
  if(!getcwd(directory,sizeof(directory))||strlen(directory)+32>=sizeof(directory))bt_fail("launch path too long");
  strcat(directory,"/.openmx-bands-XXXXXX");if(!mkdtemp(directory))bt_fail("cannot create band work directory");
  printf("CPU_TARGET_BANDS backend=SLEPc points=%d workers=%d requested=%d dimension=%d reference_Ha=%.12e full_spectrum=0 sparse_native_H_S=1 directory=%s\n",points,workers,count,dimension,target,directory);fflush(stdout);
 }
 MPI_Bcast(directory,sizeof(directory),MPI_CHAR,0,mpi_comm_level1);
 bt_parts(directory,rank,n,spin,h,im,s);
 if(!rank){
  char path[PATH_MAX];snprintf(path,sizeof(path),"%s/kpoints.bin",directory);FILE *fp=fopen(path,"wb");if(!fp)bt_fail("cannot write k points");
  int32_t nk=points;bt_write(&nk,4,1,fp);
  for(int p=1;p<=nkpath;p++)for(int j=0;j<n_perk[p];j++){
   double k[3];for(int d=1;d<=3;d++)k[d-1]=kpath[p][1][d]+(kpath[p][2][d]-kpath[p][1][d])*(double)j/(n_perk[p]-1);bt_write(k,8,3,fp);
  }
  if(fclose(fp))bt_fail("cannot close kpoint file");
 }
 MPI_Barrier(mpi_comm_level1);double snapshot_time=MPI_Wtime()-started;
 if(rank<workers){
  char helper_path[PATH_MAX],r[32],w[32],p[32],c[32],t[64],threads[32],log[PATH_MAX];
  const char *helper=cpu_band_helper(helper_path,sizeof(helper_path));if(!helper||!*helper)bt_fail("band backend is not installed");
  snprintf(r,sizeof(r),"%d",rank);snprintf(w,sizeof(w),"%d",workers);snprintf(p,sizeof(p),"%d",np);snprintf(c,sizeof(c),"%d",count);snprintf(t,sizeof(t),"%.17g",target);snprintf(threads,sizeof(threads),"%d",openmp_threads_num);
  snprintf(log,sizeof(log),"%s/worker-%06d.log",directory,rank);
  char *args[]={CpuBandOptions.python,(char*)helper,"--directory",directory,"--rank",r,"--workers",w,"--parts",p,"--count",c,"--target-ha",t,"--threads",threads,NULL};
  if(cpu_band_spawn(args,log)){fprintf(stderr,"Inspect %s\n",log);bt_fail("SLEPc worker failed; no full-spectrum fallback was performed");}
 }
 size_t total=(size_t)points*channels*count;double *local=bt_mem(total,sizeof(double)),*all=rank?NULL:bt_mem(total,sizeof(double));
 if(rank<workers){
  char path[PATH_MAX];snprintf(path,sizeof(path),"%s/values-%06d.bin",directory,rank);FILE *fp=fopen(path,"rb");if(!fp)bt_fail("missing worker energies");
  for(int point=rank;point<points;point+=workers)for(int channel=0;channel<channels;channel++){
   int32_t header[2];bt_read(header,4,2,fp);if(header[0]!=point||header[1]!=channel)bt_fail("worker point/channel ordering mismatch");
   double *values=local+((size_t)point*channels+channel)*count;bt_read(values,8,count,fp);
   for(int b=0;b<count;b++)if(!isfinite(values[b])||(b&&values[b]<values[b-1]))bt_fail("invalid returned band energies");
  }
  if(fgetc(fp)!=EOF)bt_fail("unexpected trailing worker result");fclose(fp);
 }
 MPI_Reduce(local,all,(int)total,MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);free(local);
 if(!rank){
  bt_write_band(nkpath,n_perk,kpath,kname,channels,count,points,all);
  char path[YOUSO10];snprintf(path,sizeof(path),"%s%s.Band.info",filepath,filename);FILE *fp=fopen(path,"w");if(!fp)bt_fail("cannot write band provenance");
  fprintf(fp,"Backend: SLEPc KrylovSchur shift-invert\nRequested per k/channel: %d\nMatrix dimension: %d\nReference (absolute Ha): %.17g\nReference (absolute eV): %.17g\nSelection: nearest reference, energies sorted ascending\nEigenvalues and chemical potential in .Band: Hartree\nFull spectrum: no\nOverlap: native matching S, no diagonal shift\nIndependent k workers: %d\nThreads per worker: %d\nDetailed snapshots and eigenpair audits: %s\nSnapshot wall seconds: %.9f\nBand stage wall seconds: %.9f\n",count,dimension,target,target*27.211386245988,workers,openmp_threads_num,directory,snapshot_time,MPI_Wtime()-started);
  if(fclose(fp))bt_fail("cannot close band provenance");
  printf("CPU_TARGET_BANDS_DONE requested=%d points=%d snapshot_seconds=%.6f total_seconds=%.6f normal=1\n",count,points,snapshot_time,MPI_Wtime()-started);fflush(stdout);free(all);
 }
 MPI_Barrier(mpi_comm_level1);
}

/* Keep SOC's raw H/S intact until all force consumers have finished. */
static struct {int active,nkpath,spin;int *counts;double ***path;char ***names;double *****h,*****im;double ****s;} bt_pending;
void Band_DFT_kpath(int nkpath,int *n_perk,double ***kpath,char ***kname,int spin,double *****h,double *****im,double ****s){
 if(CpuBandOptions.count&&spin==3&&CpuSOCNative_Enabled()){
  if(bt_pending.active)bt_fail("multiple projected SOC band stages are not supported");
  bt_pending.active=1;bt_pending.nkpath=nkpath;bt_pending.spin=spin;bt_pending.counts=n_perk;bt_pending.path=kpath;bt_pending.names=kname;bt_pending.h=h;bt_pending.im=im;bt_pending.s=s;
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);if(!rank){printf("CPU_TARGET_BANDS_DEFERRED until_after_force=1 operator=matching_projected_H_S\n");fflush(stdout);}return;
 }
 bt_run_bands(nkpath,n_perk,kpath,kname,spin,h,im,s);
}
void CpuTargetBands_RunPending(void){
 if(!bt_pending.active)return;
 CpuSOCNative_PrepareBands();bt_pending.active=0;
 bt_run_bands(bt_pending.nkpath,bt_pending.counts,bt_pending.path,bt_pending.names,bt_pending.spin,bt_pending.h,bt_pending.im,bt_pending.s);
}
