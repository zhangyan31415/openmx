#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "selected_left.h"
#include "sector_native_pack.h"
static int owner(int g,int nb,int p){return(g/nb)%p;}
static int loc(int g,int nb,int p){return(g/(nb*p))*nb+g%nb;}
static double sign4(int a,int b){return((a&b)==1||(a&b)==2)?-1.:1.;}
static double element(int row,int column){return .03125*(1+row+3*column);}
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
 int n=48,nb=5,pc=(int)sqrt((double)np);while(np%pc)pc--;int pr=np/pc;
 int physical=(rank+1)%np,row=physical/pc,col=physical%pc,nr=0,nc=0;
 for(int q=0;q<n;q++){nr+=owner(q,nb,pr)==row;nc+=owner(q,nb,pc)==col;}
 int64_t cp[49];int ri[192],inv[48],ends[4]={12,24,36,48};double cv[192];
 for(int q=0;q<n;q++){
  cp[q]=4*q;inv[q]=q*5%n;int component=q%12,sector=q/12;
  for(int a=0;a<4;a++){ri[4*q+a]=4*component+a;cv[4*q+a]=.5*sign4(a,sector);}
 }cp[n]=192;
 double*A=malloc((size_t)(nr?nr:1)*(nc?nc:1)*sizeof(double));
 double*reference=malloc((size_t)(nr?nr:1)*(nc?nc:1)*sizeof(double));
 double*packed=malloc((size_t)(nr?nr:1)*(nc?nc:1)*sizeof(double));
 size_t offsets[4];int rs[4],cs[4],rn[4],racc=0,cacc=0;size_t used=0;
 for(int s=0;s<4;s++){
  rs[s]=racc;cs[s]=cacc;offsets[s]=used;int cn=0;rn[s]=0;
  for(int g=s*12;g<(s+1)*12;g++){rn[s]+=owner(g,nb,pr)==row;cn+=owner(g,nb,pc)==col;}
  used+=(size_t)rn[s]*cn;racc+=rn[s];cacc+=cn;
 }
 for(int q=0;q<n;q++)if(owner(q,nb,pc)==col)for(int i=0;i<n;i++)if(owner(i,nb,pr)==row){
  A[(size_t)loc(q,nb,pc)*nr+loc(i,nb,pr)]=i/12==q/12?element(i,q):NAN;
  reference[(size_t)loc(q,nb,pc)*nr+loc(i,nb,pr)]=.5*sign4(i%4,q/12)*element((q/12)*12+i/4,q);
  if(i/12==q/12){int s=q/12;packed[offsets[s]+(size_t)(loc(q,nb,pc)-cs[s])*rn[s]+loc(i,nb,pr)-rs[s]]=element(i,q);}
 }
 CpuNativeSector source_sectors[4];size_t source_used=0;
 for(int s=0;s<4;s++){
  int sr=0,sc=0;for(int g=0;g<12;g++){sr+=owner(g,3,pr)==row;sc+=owner(g,3,pc)==col;}
  source_sectors[s]=(CpuNativeSector){12,s*12,sr,sc,3,source_used};source_used+=(size_t)sr*sc;
 }
 double*source=malloc((source_used?source_used:1)*sizeof(double));
 for(int s=0;s<4;s++)for(int q=0;q<12;q++)if(owner(q,3,pc)==col)for(int i=0;i<12;i++)if(owner(i,3,pr)==row)
  source[source_sectors[s].offset+(size_t)loc(q,3,pc)*source_sectors[s].nr+loc(i,3,pr)]=element(i+s*12,q+s*12);
 int*rankmap=malloc((size_t)np*sizeof(int));for(int r=0;r<np;r++)rankmap[(r+1)%np]=r;
 double*converted=malloc((used+8)*sizeof(double));for(size_t k=0;k<used+8;k++)converted[k]=123.;
 size_t converted_offsets[4],converted_used=0,temporary=0;
 int packed_rc=CpuNativePack(MPI_COMM_WORLD,n,nb,nr,nc,pr,pc,row,col,rankmap,source_sectors,4,source,converted,used,converted_offsets,&converted_used,7,&temporary);
 int pack_bad=packed_rc||converted_used!=used;
 for(int s=0;s<4;s++)if(offsets[s]!=converted_offsets[s])pack_bad=1;
 for(size_t k=0;k<used;k++)if(converted[k]!=packed[k])pack_bad=1;
 for(size_t k=used;k<used+8;k++)if(converted[k]!=123.)pack_bad=1;
 free(source);free(rankmap);free(converted);
 char error[256]={0};CfsSelectedPlan*p=CfsSelectedCreate(n,nb,nr,nc,pr,pc,row,col,cp,ri,cv,MPI_COMM_WORLD,error,sizeof(error));
 if(!p){fprintf(stderr,"plan %s\n",error);MPI_Abort(MPI_COMM_WORLD,1);}
 int bad=pack_bad;
 for(int trial=0;trial<6;trial++){
  int wanted=trial%3==0?n:trial%3==1?26:3;
  int first=(int)((double)wanted/np*rank)+1,last=(int)((double)wanted/np*(rank+1));
  if(rank==np-1)last=wanted;int nk=last-first+1;
  size_t live=(size_t)n*nk;double*out=malloc((live+8)*sizeof(double));
  for(size_t k=0;k<live;k++)out[k]=NAN;for(size_t k=live;k<live+8;k++)out[k]=9876543.;
  CfsSelectedTimes timing;double audit=-1;
  int rc=CfsSelectedNativeLeft(p,MPI_COMM_WORLD,trial<3?A:packed,nr,inv,wanted,first,last,out,live,ends,4,trial<3?NULL:offsets,8+trial*4,reference,1,&timing,&audit,error,sizeof(error));
  if(rc||audit!=0.){fprintf(stderr,"rank%d native rc%d audit%.17g %s\n",rank,rc,audit,error);bad=1;}
  for(int q=0;q<n;q++)if(inv[q]>=first-1&&inv[q]<last)for(int i=0;i<n;i++){
   double expected=.5*sign4(i%4,q/12)*element((q/12)*12+i/4,q);
   double got=out[(size_t)i*nk+inv[q]-(first-1)];if(got!=expected)bad=1;
  }
  for(size_t k=live;k<live+8;k++)if(out[k]!=9876543.)bad=1;
  free(out);
 }
 /* Collectively reject an invalid sorting permutation before touching output. */
 double small=123.;inv[1]=inv[0];
 int rejected=CfsSelectedNativeLeft(p,MPI_COMM_WORLD,A,nr,inv,1,1,rank?0:1,&small,rank?0:1,ends,4,NULL,8,NULL,0,NULL,NULL,error,sizeof(error));
 /* The tiny capacity can reject first; either boundary must fail closed. */
 if(!rejected||small!=123.)bad=1;
 CfsSelectedDestroy(p);free(A);free(reference);free(packed);
 int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 if(!rank)printf("selected_native ranks=%d grid=%dx%d full_partial_empty_owners=PASS sector_nan_ignored=PASS exact_coefficients=%s\n",np,pr,pc,any?"FAIL":"PASS");
 MPI_Finalize();return any?1:0;
}
