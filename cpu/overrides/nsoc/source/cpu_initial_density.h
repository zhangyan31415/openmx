#ifndef CPU_INITIAL_DENSITY_H
#define CPU_INITIAL_DENSITY_H
#include <limits.h>
#include "cpu_initial_grid_model.h"

static void CpuInitialDensityFail(const char *why)
{fprintf(stderr,"CPU_INITIAL_DENSITY rejects %s\n",why);MPI_Abort(mpi_comm_level1,98);exit(98);}

static void CpuInitialDensityBroadcast(double *p,size_t count)
{
  while(count){int n=count>(size_t)INT_MAX?INT_MAX:(int)count;
    MPI_Bcast(p,n,MPI_DOUBLE,Host_ID,mpi_comm_level1);p+=n;count-=n;
  }
}

static void CpuInitialDensityApply(void)
{
  static int done=0;
  const char *path=getenv("OPENMX_CPU_NSOC_INITIAL_GRID");
  int enabled=path&&*path&&strcmp(path,"0"),lo,hi,rank,np;
  MPI_Comm_rank(mpi_comm_level1,&rank);MPI_Comm_size(mpi_comm_level1,&np);
  MPI_Allreduce(&enabled,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&enabled,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo!=hi)CpuInitialDensityFail("inconsistent selection across MPI ranks");
  if(!enabled||done)return;
  if(SpinP_switch!=0||SO_switch||MD_switch||Solver!=2||Cnt_switch||LNO_flag||
     Scf_RestartFromFile!=0||SCF_Control_Temp)
    CpuInitialDensityFail("fixed fresh scalar Cluster NSOC required");
  unsigned long long hash=1469598103934665603ULL,hashlo,hashhi;
  for(const unsigned char *p=(const unsigned char*)path;*p;p++)hash=(hash^*p)*1099511628211ULL;
  MPI_Allreduce(&hash,&hashlo,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&hash,&hashhi,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,mpi_comm_level1);
  if(hashlo!=hashhi)CpuInitialDensityFail("inconsistent initial-grid path");
  double start=MPI_Wtime();int normalize=!strcmp(path,"normalize");
  CpuInitialGridModel model;memset(&model,0,sizeof(model));
  if(!normalize){
    int status=0,dims[6];
    if(rank==Host_ID)status=CpuInitialGridRead(path,&model);
    MPI_Bcast(&status,1,MPI_INT,Host_ID,mpi_comm_level1);
    if(status)CpuInitialDensityFail("invalid template file");
    if(rank==Host_ID){for(int i=0;i<3;i++)dims[i]=model.n[i];dims[3]=model.frame_count;dims[4]=model.warp_count;dims[5]=model.isolated_z;}
    MPI_Bcast(dims,6,MPI_INT,Host_ID,mpi_comm_level1);
    size_t count=(size_t)dims[0]*dims[1]*dims[2],frames=(size_t)dims[3]*4;
    int bad=0,anybad;
    if(rank!=Host_ID){
      for(int i=0;i<3;i++)model.n[i]=dims[i];model.frame_count=dims[3];
      model.frames=malloc(frames*sizeof(double));model.values=malloc(count*sizeof(double));
      bad=!model.frames||!model.values;
      model.warp_count=dims[4];
      model.isolated_z=dims[5];
      if(model.warp_count){
        model.wave=malloc((size_t)model.warp_count*2*sizeof(int));
        model.warp_coeff=malloc((size_t)model.frame_count*model.warp_count*6*sizeof(double));
        bad|=!model.wave||!model.warp_coeff;
      }
    }
    MPI_Allreduce(&bad,&anybad,1,MPI_INT,MPI_MAX,mpi_comm_level1);
    if(anybad)CpuInitialDensityFail("template workspace allocation");
    MPI_Bcast(model.origin,3,MPI_DOUBLE,Host_ID,mpi_comm_level1);
    MPI_Bcast(model.inverse,9,MPI_DOUBLE,Host_ID,mpi_comm_level1);
    CpuInitialDensityBroadcast(model.frames,frames);CpuInitialDensityBroadcast(model.values,count);
    if(model.warp_count){
      MPI_Bcast(model.warp_inverse,4,MPI_DOUBLE,Host_ID,mpi_comm_level1);
      MPI_Bcast(model.wave,2*model.warp_count,MPI_INT,Host_ID,mpi_comm_level1);
      CpuInitialDensityBroadcast(model.warp_coeff,(size_t)model.frame_count*model.warp_count*6);
    }
  }
  unsigned long long n2d=(unsigned long long)Ngrid1*Ngrid2;
  unsigned long long begin=((rank*n2d+np-1)/np)*(unsigned long long)Ngrid3;
  double local[4]={0.,0.,0.,0.},global[4],expected=-system_charge;
  int bad=0,anybad;
  for(int a=1;a<=atomnum;a++)expected+=Spe_Core_Charge[WhatSpecies[a]];
  for(int b=0;b<My_NumGridB_AB;b++){
    double rho=Density_Grid_B[0][b]+Density_Grid_B[1][b],delta=0.;local[0]+=rho;
    if(!normalize){
      unsigned long long gn=begin+(unsigned long long)b;
      unsigned long long ij=gn/(unsigned long long)Ngrid3;
      unsigned long long i=ij/(unsigned long long)Ngrid2,j=ij%(unsigned long long)Ngrid2;
      unsigned long long k=gn%(unsigned long long)Ngrid3;double r[3];
      for(int a=0;a<3;a++)r[a]=Grid_Origin[a+1]+i*gtv[1][a+1]+j*gtv[2][a+1]+k*gtv[3][a+1];
      delta=CpuInitialGridSample(&model,r);
    }
    local[2]+=delta;rho+=delta;if(!isfinite(rho)){bad=1;rho=0.;}
    if(rho<0.){local[3]-=rho;rho=0.;}
    local[1]+=rho;Density_Grid_B[0][b]=Density_Grid_B[1][b]=.5*rho;
  }
  MPI_Allreduce(&bad,&anybad,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
  double integral=global[1]*GridVol;
  if(anybad||!isfinite(integral)||integral<=0.||!isfinite(expected)||expected<=0.)
    CpuInitialDensityFail("nonfinite or empty initial electron density");
  double scale=expected/integral;
  for(int b=0;b<My_NumGridB_AB;b++){
    Density_Grid_B[0][b]*=scale;Density_Grid_B[1][b]*=scale;
  }
  if(rank==Host_ID){printf("CPU_INITIAL_DENSITY mode=%s reference_atomic_density_unchanged=1 original_electrons=%.12e correction_electrons=%.12e clipped_electrons=%.12e scale=%.12e final_electrons=%.12e seconds=%.6f\n",normalize?"normalize":"layer_bond_correction",global[0]*GridVol,global[2]*GridVol,global[3]*GridVol,scale,expected,MPI_Wtime()-start);fflush(stdout);}
  CpuInitialGridFree(&model);done=1;
}
#endif
