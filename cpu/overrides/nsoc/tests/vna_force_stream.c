/* Compare the streamed scalar/spinor-trace force to unmodified native dHVNA.
 * Build with the extracted native function in dHVNA_reference.inc. */
#include <stdio.h>
#include "openmx_common.h"
#include "cpu_vna_shared.h"
#include "cpu_vna_stream.h"
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
static int p=4;
int flag_export_DM;
int CpuGridShared_Requested(void){return 1;}
double dampingF(double cut,double r){(void)cut;return r>0.?.75:1.;}
double deri_dampingF(double cut,double r){(void)cut;(void)r;return .25;}
static float value(int a,int q,int m,int l){return (float)(16*a+4*q+2*m+l)/64.f;}
static float gradient(int a,int q,int m,int l,int xyz){return a==q?0.f:(float)((xyz+1)*(3*a+5*q+2*m+l))/128.f;}
static double **matrix(int n,int m){double **v=calloc(n,sizeof(*v));for(int i=0;i<n;i++)v[i]=calloc(m,sizeof(**v));return v;}
static double ****blocks(int locals,int neighbors){double ****v=calloc(locals,sizeof(*v));for(int a=0;a<locals;a++){v[a]=calloc(neighbors,sizeof(**v));for(int h=0;h<neighbors;h++)v[a][h]=matrix(2,2);}return v;}
static double *****components(int count,int locals,int neighbors){double *****v=calloc(count,sizeof(*v));for(int k=0;k<count;k++)v[k]=blocks(locals,neighbors);return v;}
/* Mock only the expensive native integral generator; the production force
 * consumer and the reference dHVNA are compiled without alteration. */
double CpuVnaDeferred_Stream(float *****view,CpuVnaPairGradient callback,void *opaque){
 (void)view;int a=M2G[1];if(CpuVnaShared_Raw0Dropped()){
  for(int h=0;h<atomnum;h++)for(int m=0;m<2;m++)for(int l=0;l<p;l++)CpuVnaShared_Row(0,a,h,m)[l]=value(a,natn[a][h],m,l);CpuVnaShared_RestoreRaw0();
 }
 float g[3*2*4];for(int h=0;h<atomnum;h++){
  int q=natn[a][h];for(int xyz=0;xyz<3;xyz++)for(int m=0;m<2;m++)for(int l=0;l<p;l++)g[xyz*2*p+m*p+l]=gradient(a,q,m,l,xyz);callback(opaque,1,h,g,2,p);
 }return 0.;
}
#include "dHVNA_reference.inc"
int main(int argc,char **argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);assert(np>=2);omp_set_num_threads(1);
 int soc=getenv("OPENMX_TEST_SOC")!=NULL;
 mpi_comm_level1=MPI_COMM_WORLD;atomnum=np;Matomnum=1;MatomnumF=np-1;MD_switch=0;Solver=2;SpinP_switch=soc?3:0;SO_switch=soc?1:0;DFTSCF_loop=1;Cnt_switch=0;ProExpn_VNA=1;F_VNA_flag=1;scf_stress_flag=0;List_YOUSO[34]=1;List_YOUSO[35]=1;
 if(soc)setenv("OPENMX_CPU_SOC_MEMORY","1",1);
 Mixing_switch=1;Scf_RestartFromFile=soc?1:0;
 int *owners=calloc(np+1,sizeof(int)),*norb=calloc(np+1,sizeof(int));FNAN=calloc(np+1,sizeof(int));natn=calloc(np+1,sizeof(int*));ncn=calloc(np+1,sizeof(int*));WhatSpecies=calloc(np+1,sizeof(int));Gxyz=matrix(np+1,64);Dis=matrix(np+1,np);atv=matrix(1,4);Spe_Total_CNO=calloc(2,sizeof(int));Spe_Total_NO=calloc(2,sizeof(int));Spe_Atom_Cut1=calloc(2,sizeof(double));
 for(int s=0;s<2;s++){Spe_Total_CNO[s]=Spe_Total_NO[s]=2;Spe_Atom_Cut1[s]=10.;}
 for(int a=1;a<=np;a++){owners[a]=a-1;norb[a]=2;FNAN[a]=np-1;natn[a]=calloc(np,sizeof(int));ncn[a]=calloc(np,sizeof(int));WhatSpecies[a]=a%2;Gxyz[a][1]=a*.25;Gxyz[a][2]=a*a*.125;for(int h=0;h<np;h++)natn[a][h]=(a+h-1)%np+1;}
 for(int a=1;a<=np;a++)for(int h=0;h<np;h++){int b=natn[a][h];Dis[a][h]=hypot(Gxyz[a][1]-Gxyz[b][1],Gxyz[a][2]-Gxyz[b][2]);}
 M2G=calloc(2,sizeof(int));M2G[1]=rank+1;F_M2G=calloc(np+1,sizeof(int));F_G2M=calloc(np+1,sizeof(int));F_M2G[1]=rank+1;F_G2M[rank+1]=1;int at=2;for(int a=1;a<=np;a++)if(a!=rank+1){F_M2G[at]=a;F_G2M[a]=at++;}
 RMI1=calloc(2,sizeof(int**));RMI1[1]=calloc(np,sizeof(int*));for(int h=0;h<np;h++){RMI1[1][h]=calloc(np,sizeof(int));for(int q=0;q<np;q++)RMI1[1][h][q]=(natn[rank+1][q]-natn[rank+1][h]+np)%np;}
 VNA_proj_ene=calloc(2,sizeof(double**));for(int s=0;s<2;s++){VNA_proj_ene[s]=matrix(2,1);VNA_proj_ene[s][0][0]=.5*(s+1);VNA_proj_ene[s][1][0]=.25*(s+1);}
 HVNA=blocks(2,np);HVNA2=components(4,2,np);HVNA3=components(4,2,np);double *****density=components(soc?4:1,np+1,np);
 for(int local=1;local<=np;local++){int a=F_M2G[local];for(int h=0;h<np;h++){int b=natn[a][h];for(int i=0;i<2;i++)for(int j=0;j<2;j++){
  density[0][local][h][i][j]=(double)(7*(a+b)+3*(i+j))/128.;
  if(soc){density[1][local][h][i][j]=(double)(5*(a+b)+2*(i+j)+1)/256.;density[2][local][h][i][j]=.125;density[3][local][h][i][j]=-.25;}
 }}}
 for(int h=0;h<np;h++)for(int i=0;i<2;i++)for(int j=0;j<2;j++){HVNA[1][h][i][j]=(double)(3*h+2*i+j+1)/64.;for(int k=1;k<4;k++){HVNA2[k][1][h][i][j]=(double)(k+h+i+j)/512.;HVNA3[k][1][h][i][j]=(double)(2*k+h+i+j)/512.;}}
 int image0[4]={0},*images[1]={image0};setenv("OPENMX_CPU_VNA_SHARED","1",1);setenv("OPENMX_CPU_VNA_FORCE_STREAM",argc>1?"2":"1",1);
 CpuVnaShared_Init(MPI_COMM_WORLD,np,owners,FNAN,natn,ncn,images,norb,p,np-1,2);CpuVnaShared_ConfigureDeferred(1);int a=rank+1;
 for(int h=0;h<np;h++)for(int m=0;m<2;m++)for(int l=0;l<p;l++)CpuVnaShared_Row(0,a,h,m)[l]=value(a,natn[a][h],m,l);CpuVnaShared_Publish(0);if(argc>1)CpuVnaShared_DropRaw0();
 DS_VNA=calloc(4,sizeof(*DS_VNA));for(int k=0;k<4;k++){DS_VNA[k]=calloc(3,sizeof(**DS_VNA));for(int local=1;local<=2;local++){DS_VNA[k][local]=calloc(np,sizeof(***DS_VNA));for(int h=0;h<np;h++){DS_VNA[k][local][h]=calloc(2,sizeof(****DS_VNA));for(int m=0;m<2;m++)if(local==1)DS_VNA[k][local][h][m]=CpuVnaShared_Row(k,a,h,m);}}}
 assert(CpuVnaForceStream_Run(density));assert(CpuVnaShared_Raw0Dropped());assert(CpuVnaShared_Active());double got[3];for(int xyz=0;xyz<3;xyz++)got[xyz]=Gxyz[a][41+xyz];
 /* Production has no more readers after the stream force. Rebuild values
  * only to evaluate the independent native-reference path below. */
 for(int h=0;h<np;h++)for(int m=0;m<2;m++)for(int l=0;l<p;l++)CpuVnaShared_Row(0,a,h,m)[l]=value(a,natn[a][h],m,l);
 CpuVnaShared_RestoreRaw0();assert(!CpuVnaShared_Raw0Dropped());
 /* Populate the complete weighted cache only for native-reference evaluation. */
 for(int h=0;h<np;h++)for(int m=0;m<2;m++)for(int l=0;l<p;l++)for(int xyz=0;xyz<3;xyz++){int q=natn[a][h],angular=l==0?0:1;CpuVnaShared_Row(xyz+1,a,h,m)[l]=(float)(VNA_proj_ene[WhatSpecies[q]][angular][0]*gradient(a,q,m,l,xyz));}
 CpuVnaShared_DerivativesReady();CpuVnaShared_BeginWeighted();CpuVnaShared_Publish(1);double expected[3]={0.,0.,0.};double **x=matrix(2,2),**y=matrix(2,2),**z=matrix(2,2),**axis[3]={x,y,z};
 for(int q=0;q<np;q++){
  CpuVnaShared_BindAtom(DS_VNA,2,natn[a][q]);dHVNA(0,1,0,q,DS_VNA,HVNA2,HVNA3,x,y,z);double pref=q?4.:2.;for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int xyz=0;xyz<3;xyz++){double d=soc?.5*(density[0][1][q][i][j]+density[1][1][q][i][j]):density[0][1][q][i][j];expected[xyz]+=pref*d*axis[xyz][i][j];}
 }
 CpuVnaShared_BindCenter(DS_VNA,2,a,RMI1[1]);for(int h=0;h<np;h++)for(int q=h;q<np;q++){
  dHVNA(1,1,h,q,DS_VNA,HVNA2,HVNA3,x,y,z);double pref=q==h?2.:4.;int row=F_G2M[natn[a][h]],slot=RMI1[1][h][q];for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int xyz=0;xyz<3;xyz++){double d=soc?.5*(density[0][row][slot][i][j]+density[1][row][slot][i][j]):density[0][row][slot][i][j];expected[xyz]+=pref*d*axis[xyz][i][j];}
 }
 double error=0.;for(int xyz=0;xyz<3;xyz++){double e=fabs(expected[xyz]-got[xyz]);if(e>error)error=e;}double maximum;MPI_Allreduce(&error,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);printf("rank=%d native=(%.15g %.15g %.15g) stream=(%.15g %.15g %.15g) error=%.3g\n",rank,expected[0],expected[1],expected[2],got[0],got[1],got[2],error);assert(maximum<1.e-11);CpuVnaShared_Destroy();if(!rank)printf("PASS streamed force native reference maxerr=%.3g\n",maximum);MPI_Finalize();return 0;
}
