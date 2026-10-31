#ifndef GSOC_NATIVE_DENSITY_AUDIT_H
#define GSOC_NATIVE_DENSITY_AUDIT_H
#include <stdint.h>
#include <math.h>
/* Tiny integration oracle only. No collective or extra dense matrix, and
 * disabled in every application timing case. Each rank writes its own blocks. */
static void gsoc_native_density_audit(const double *ko)
{
  const char *flag=getenv("OPENMX_CPU_GSOC_DENSITY_AUDIT");
  if(!flag||!strcmp(flag,"0"))return;
  int rank,np,N=0;
  MPI_Comm_rank(mpi_comm_level1,&rank);MPI_Comm_size(mpi_comm_level1,&np);
  for(int a=1;a<=atomnum;a++)N+=Spe_Total_CNO[WhatSpecies[a]];
  if(strcmp(flag,"1")||atomnum>32||N>256){fprintf(stderr,"GSOC_DENSITY_AUDIT restricted to tiny integration tests\n");MPI_Abort(mpi_comm_level1,96);}
  char name[96];snprintf(name,sizeof(name),"gsoc_density_rank%d.bin",rank);
  FILE *file=fopen(name,"wb");
  if(!file){fprintf(stderr,"GSOC_DENSITY_AUDIT open failed\n");MPI_Abort(mpi_comm_level1,96);}
  int32_t header[8]={0x47534441,1,rank,np,Matomnum,atomnum,N,2*N};
  if(fwrite(header,sizeof(int32_t),8,file)!=8)MPI_Abort(mpi_comm_level1,96);
  double therm[2]={ChemP,Beta};fwrite(therm,sizeof(double),2,file);
  for(int band=1;band<=2*N;band++){
    double x=(ko[band]-ChemP)*Beta;
    if(x < -60.)x=-60.;if(x>60.)x=60.;
    double f=FermiFunc_NC(x,band),pair[2]={ko[band],f>1e-13?f:0.};
    fwrite(pair,sizeof(double),2,file);
  }
  double traces[2]={0.,0.};
  for(int a=1;a<=Matomnum;a++){
    int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];
    int32_t atom[3]={ga,na,FNAN[ga]+1};fwrite(atom,sizeof(int32_t),3,file);
    for(int h=0;h<=FNAN[ga];h++){
      int gb=natn[ga][h],nb=Spe_Total_CNO[WhatSpecies[gb]];
      int32_t block[2]={gb,nb};fwrite(block,sizeof(int32_t),2,file);
      for(int i=0;i<na;i++)for(int j=0;j<nb;j++){
        double v[11];for(int s=0;s<4;s++)v[s]=DM[0][s][a][h][i][j];
        for(int s=0;s<2;s++)v[4+s]=iDM[0][s][a][h][i][j];
        for(int s=0;s<4;s++)v[6+s]=EDM[s][a][h][i][j];
        v[10]=OLP[0][a][h][i][j];
        for(int s=0;s<11;s++)if(!isfinite(v[s]))MPI_Abort(mpi_comm_level1,96);
        traces[0]+=(v[0]+v[1])*v[10];traces[1]+=(v[6]+v[7])*v[10];
        if(fwrite(v,sizeof(double),11,file)!=11)MPI_Abort(mpi_comm_level1,96);
      }
    }
  }
  fwrite(traces,sizeof(double),2,file);
  if(fclose(file))MPI_Abort(mpi_comm_level1,96);
}
#endif
