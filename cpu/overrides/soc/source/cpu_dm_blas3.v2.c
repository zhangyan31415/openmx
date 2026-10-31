#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mkl.h>
#include "cpu_dm_blas3.v2.h"

int cpu_dm_blocks(int n,int bands,const void *vectors,const double *w,
                  int count,const int *blocks,int entries,int edm,int tile,double *out)
{
  const MKL_Complex16 *v=(const MKL_Complex16*)vectors;
  int p,c,first,last,m,nc,maxm=0,maxc=0,i,j,k,t,nb,pair,base,b,offset;
  int components=edm?4:6;
  MKL_Complex16 *A=NULL,*B=NULL,*G=NULL;
  const MKL_Complex16 one={1.0,0.0};
  if(n<0||bands<0||count<0||entries<0||tile<=0||!(edm==0||edm==1))return -1;
  if((count&&!blocks)||(entries&&!out)||(bands&&(!v||!w)))return -1;
  for(c=0;c<count;c++){
    const int *x=blocks+5*c;
    if(x[0]<0||x[1]<=0||x[2]<0||x[3]<=0||x[4]<0||
       (int64_t)x[0]+x[1]>n||(int64_t)x[2]+x[3]>n||
       (int64_t)x[4]+(int64_t)x[1]*x[3]>entries)return -1;
  }
  if(entries)memset(out,0,(size_t)components*entries*sizeof(double));
  if(!bands||!count)return 0;
  for(k=0;k<bands&&w[k]==0.0;k++);
  if(k==bands)return 0;
  for(first=0;first<count;first=last){
    m=blocks[5*first+1];nc=0;last=first;
    while(last<count&&blocks[5*last]==blocks[5*first]&&blocks[5*last+1]==m){
      if((int64_t)nc+blocks[5*last+3]>2147483647)return -1;
      nc+=blocks[5*last+3];last++;
    }
    if(m>maxm)maxm=m;if(nc>maxc)maxc=nc;
  }
  if(tile>bands)tile=bands;
  A=(MKL_Complex16*)malloc((size_t)maxm*tile*sizeof(*A));
  B=(MKL_Complex16*)malloc((size_t)maxc*tile*sizeof(*B));
  G=(MKL_Complex16*)malloc((size_t)maxm*maxc*sizeof(*G));
  if(!A||!B||!G){free(A);free(B);free(G);return -2;}
  for(first=0;first<count;first=last){
    m=blocks[5*first+1];nc=0;last=first;
    while(last<count&&blocks[5*last]==blocks[5*first]&&blocks[5*last+1]==m){nc+=blocks[5*last+3];last++;}
    for(pair=0;pair<3;pair++){
      int sa=pair==1?n:0,sb=pair==0?0:n;
      memset(G,0,(size_t)m*nc*sizeof(*G));
      for(k=0;k<bands;k+=tile){
        nb=bands-k<tile?bands-k:tile;
        for(t=0;t<nb;t++){
          for(i=0;i<m;i++){
            MKL_Complex16 z={0,0};
            if(w[k+t]!=0){z=v[(size_t)(k+t)*2*n+sa+blocks[5*first]+i];z.real*=w[k+t];z.imag*=w[k+t];}
            A[(size_t)t*m+i]=z;
          }
          base=0;
          for(b=first;b<last;b++){
            for(j=0;j<blocks[5*b+3];j++){
              MKL_Complex16 z={0,0};
              if(w[k+t]!=0)z=v[(size_t)(k+t)*2*n+sb+blocks[5*b+2]+j];
              B[(size_t)t*nc+base+j]=z;
            }
            base+=blocks[5*b+3];
          }
        }
        cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,m,nc,nb,&one,A,m,B,nc,&one,G,m);
      }
      base=0;
      for(b=first;b<last;b++){
        int q=blocks[5*b+3];offset=blocks[5*b+4];
        for(i=0;i<m;i++)for(j=0;j<q;j++){
          MKL_Complex16 z=G[i+(size_t)(base+j)*m];p=offset+i*q+j;
          if(pair==0){out[p]=z.real;if(!edm)out[(size_t)4*entries+p]=-z.imag;}
          if(pair==1){out[(size_t)entries+p]=z.real;if(!edm)out[(size_t)5*entries+p]=-z.imag;}
          if(pair==2){out[(size_t)2*entries+p]=z.real;out[(size_t)3*entries+p]=-z.imag;}
        }
        base+=q;
      }
    }
  }
  free(A);free(B);free(G);return 0;
}

#include <math.h>
int cpu_dm_compare(int count,const double *a,const double *b,double *max_abs,double *scale)
{
  int i;double err=0.0,mag=0.0;
  if(count<0||!max_abs||!scale||(count&&(!a||!b)))return -1;
  for(i=0;i<count;i++){
    double d;
    if(!isfinite(a[i])||!isfinite(b[i]))return -2;
    d=fabs(a[i]-b[i]);
    if(!isfinite(d))return -2;
    if(d>err)err=d;if(fabs(a[i])>mag)mag=fabs(a[i]);
  }
  *max_abs=err;*scale=mag;return 0;
}
