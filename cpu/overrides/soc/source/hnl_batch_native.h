#include "hnl_batch_kernel.h"
#define BATCH_SO_OUT(s) &(s)[0],&(s)[1],&(s)[2],&(s)[3],&(s)[4],&(s)[5],&(s)[6],&(s)[7],&(s)[8],&(s)[9],&(s)[10],&(s)[11],&(s)[12],&(s)[13],&(s)[14],&(s)[15],&(s)[16],&(s)[17]
typedef struct {double*A,*B,*C;unsigned long long calls[5],audits[5],checked;double max_error,max_reference;} CpuHnlBatchWork;
static HnlWeights CpuHnlBW;
static CpuHnlBatchWork*CpuHnlBP;
static int CpuHnlBM=0,CpuHnlBT=0,CpuHnlBN=0,CpuHnlBK=0;
static size_t CpuHnlBBytes=0;
static void CpuHnlBatchOracle(int ell,int branch,int a,int b,double out[3]){
 double v[2][4][7]={0};double *row[2][4][1];double **nei[2][4][1];double ***atom[2][4][1];double ****field[2][4];double *****spin[2];
 for(int j=0;j<2;j++){spin[j]=field[j];for(int d=0;d<4;d++){field[j][d]=atom[j][d];atom[j][d][0]=nei[j][d];nei[j][d][0]=row[j][d];row[j][d][0]=v[j][d];}}
 v[branch][1][a]=1.;v[branch][0][b]=1.;double s[18]={0};int l=0;
 dHNL_SO(BATCH_SO_OUT(s),1.,0.,0.,branch==0?1.:0.,branch==1?1.:0.,2*ell,&l,0,0,0,0,0,0,spin);
 if(l!=2*ell+1||s[0]!=0.||s[3]!=0.||s[12]!=-s[9])CpuF3TileFail("HNL_BATCH angular oracle contract");
 out[0]=s[9];out[1]=s[6];out[2]=s[15];
}
static int CpuHnlBatchEligible(int species){
 if(!CpuHnlBM||SpinP_switch!=3||Cnt_switch!=0||VPS_j_dependency[species]!=1)return 0;
 int K=0;for(int r=1;r<=Spe_Num_RVPS[species];r++){int ell=Spe_VPS_List[species][r];if(ell<0||ell>3)return 0;K+=2*ell+1;}
 return K>0&&K==Spe_Total_VPS_Pro[species]&&K<=CpuHnlBK;
}
static void CpuHnlBatchStart(void){
 const char*e=getenv("OPENMX_CPU_HNL_BATCH");int mode=!e||!strcmp(e,"0")?0:!strcmp(e,"1")?1:!strcmp(e,"2")?2:-1,in[2]={mode,-mode},out[2];
 MPI_Allreduce(in,out,2,MPI_INT,MPI_MIN,mpi_comm_level1);if(out[0]<0||out[0]!=-out[1])CpuF3TileFail("HNL_BATCH must agree0/1/2 across ranks");
 if(CpuHnlBP)CpuF3TileFail("HNL_BATCH nested lifetime");CpuHnlBM=mode;if(!mode)return;
 int okay=SpinP_switch==3&&Cnt_switch==0&&!omp_in_parallel(),all;MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)CpuF3TileFail("HNL_BATCH requires serial noncontracted SOC Force_HNL");
 CpuHnlBT=omp_get_max_threads();CpuHnlBN=List_YOUSO[7];CpuHnlBK=0;for(int s=0;s<SpeciesNum;s++)if(Spe_Total_VPS_Pro[s]>CpuHnlBK)CpuHnlBK=Spe_Total_VPS_Pro[s];
 if(CpuHnlBT<1||CpuHnlBT>256||CpuHnlBN<1||CpuHnlBN>256||CpuHnlBK<1||CpuHnlBK>4096)CpuF3TileFail("HNL_BATCH bounded dimension admission");
 CpuHnlBP=calloc(CpuHnlBT,sizeof(*CpuHnlBP));okay=CpuHnlBP!=NULL;MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)CpuF3TileFail("HNL_BATCH metadata allocation");
 size_t na=(size_t)3*CpuHnlBN*2*CpuHnlBK,nb=(size_t)2*CpuHnlBK*4*CpuHnlBN,nc=(size_t)3*CpuHnlBN*4*CpuHnlBN;
 okay=1;for(int t=0;t<CpuHnlBT;t++){CpuHnlBP[t].A=mkl_malloc(na*8,64);CpuHnlBP[t].B=mkl_malloc(nb*8,64);CpuHnlBP[t].C=mkl_malloc(nc*8,64);if(!CpuHnlBP[t].A||!CpuHnlBP[t].B||!CpuHnlBP[t].C)okay=0;}
 MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)CpuF3TileFail("HNL_BATCH panel allocation");
 CpuHnlBBytes=(na+nb+nc)*8*CpuHnlBT;hnl_weights_init(&CpuHnlBW,CpuHnlBatchOracle);
}
static void CpuHnlBatchStop(void){
 if(!CpuHnlBM)return;unsigned long long calls[5]={0},audits=0,checked=0;double error=0.,reference=0.;
 for(int t=0;t<CpuHnlBT;t++){for(int s=0;s<5;s++){calls[s]+=CpuHnlBP[t].calls[s];audits+=CpuHnlBP[t].audits[s];}checked+=CpuHnlBP[t].checked;if(CpuHnlBP[t].max_error>error)error=CpuHnlBP[t].max_error;if(CpuHnlBP[t].max_reference>reference)reference=CpuHnlBP[t].max_reference;mkl_free(CpuHnlBP[t].A);mkl_free(CpuHnlBP[t].B);mkl_free(CpuHnlBP[t].C);}
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);char line[768];int n=snprintf(line,sizeof(line),"CPU_HNL_BATCH rank=%d mode=%d threads=%d c0=%llu c1=%llu c2=%llu c3=%llu c4=%llu audits=%llu checked=%llu max_scaled_error=%.17g max_reference=%.17g panel_bytes=%zu\n",rank,CpuHnlBM,CpuHnlBT,calls[0],calls[1],calls[2],calls[3],calls[4],audits,checked,error,reference,CpuHnlBBytes);
 if(n<0||(size_t)n>=sizeof(line))CpuF3TileFail("HNL_BATCH output overflow");fwrite(line,1,n,stdout);fflush(stdout);free(CpuHnlBP);CpuHnlBP=NULL;CpuHnlBM=0;
}
static void CpuHnlBatchAdd(int site,int species,int la,int lh,int nl,int ra,int rh,int nr,double esign,double fugou,int transpose,double******DS,dcomplex***Hx,dcomplex***Hy,dcomplex***Hz){
 int t=omp_get_thread_num(),K=Spe_Total_VPS_Pro[species],radial=Spe_Num_RVPS[species];
 if(t<0||t>=CpuHnlBT||nl<1||nr<1||nl>CpuHnlBN||nr>CpuHnlBN||K>CpuHnlBK||radial<1||radial>256||site<0||site>4)CpuF3TileFail("HNL_BATCH runtime dimensions");
 CpuHnlBatchWork*p=CpuHnlBP+t;int ell[256];double ep[256],em[256];int actualK=0;for(int r=0;r<radial;r++){ell[r]=Spe_VPS_List[species][r+1];ep[r]=Spe_VNLE[0][species][r];em[r]=Spe_VNLE[1][species][r];if(ell[r]<0||ell[r]>3)CpuF3TileFail("HNL_BATCH angular admission changed");actualK+=2*ell[r]+1;}
 if(actualK!=K)CpuF3TileFail("HNL_BATCH projector column mismatch");
 double**left[2][3],**right[2];for(int b=0;b<2;b++){right[b]=DS[b][0][ra][rh];for(int d=0;d<3;d++)left[b][d]=DS[b][d+1][la][lh];}
 hnl_batch_compute(&CpuHnlBW,nl,nr,K,radial,ell,ep,em,esign,left,right,p->A,p->B,p->C);
 int audit=CpuHnlBM==2&&p->audits[site]<8&&p->calls[site]%13==0;p->calls[site]++;
 dcomplex***out[3]={Hx,Hy,Hz};
 for(int m=0;m<nl;m++)for(int n=0;n<nr;n++){
  double value[18];hnl_batch_value(p->C,nl,nr,m,n,fugou,value);
  if(audit){double ref[18]={0};int l=0;for(int r=0;r<radial;r++){int L=ell[r];dHNL_SO(BATCH_SO_OUT(ref),fugou,(L+1.)/(2*L+1.),L/(2*L+1.),esign*ep[r],esign*em[r],2*L,&l,la,lh,m,ra,rh,n,DS);}if(l!=K)CpuF3TileFail("HNL_BATCH native cursor");for(int c=0;c<18;c++){double err=fabs(value[c]-ref[c])/(1+fabs(ref[c]));if(!isfinite(err)||err>3e-12)CpuF3TileFail("HNL_BATCH native component mismatch");if(err>p->max_error)p->max_error=err;if(fabs(ref[c])>p->max_reference)p->max_reference=fabs(ref[c]);p->checked++;}}
  int i=transpose?n:m,j=transpose?m:n;
  for(int d=0;d<3;d++)for(int spin=0;spin<3;spin++){out[d][spin][i][j].r+=value[3*spin+d];out[d][spin][i][j].i+=value[9+3*spin+d];}
 }
 if(audit)p->audits[site]++;
}
#undef BATCH_SO_OUT
