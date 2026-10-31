#include <math.h>
#include <string.h>
#include "soc_spin.h"
typedef struct { double r,i; } SC;
static SC mul(SC a,SC b){SC c={a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};return c;}
static SC conjz(SC a){a.i=-a.i;return a;}
static const int rr[4]={0,1,0,1},cc[4]={0,1,1,0};
int SocSpin_Transform(const double spin[8],double transform[64]){
 SC u[2][2];int a,b,i,j,k;
 for(a=0;a<2;a++)for(b=0;b<2;b++){u[a][b].r=spin[2*(2*a+b)];u[a][b].i=spin[2*(2*a+b)+1];if(!isfinite(u[a][b].r)||!isfinite(u[a][b].i))return 1;}
 for(a=0;a<2;a++)for(b=0;b<2;b++){
  SC z={0,0};for(k=0;k<2;k++){SC t=mul(u[a][k],conjz(u[b][k]));z.r+=t.r;z.i+=t.i;}
  if(hypot(z.r-(a==b),z.i)>1e-12)return 2;
 }
 memset(transform,0,64*sizeof(double));
 for(j=0;j<8;j++)for(i=0;i<4;i++){
  SC z=mul(u[rr[i]][rr[j/2]],conjz(u[cc[i]][cc[j/2]]));
  if(j%2){double t=z.r;z.r=-z.i;z.i=t;}
  transform[8*(2*i)+j]=z.r;transform[8*(2*i+1)+j]=z.i;
 }
 return 0;
}
int SocSpin_Reverse(int field,int *source,double *sign){
 static const int reverse[8]={0,1,2,3,6,7,4,5};
 if(field<0||field>=8||!source||!sign)return 1;
 *source=reverse[field];*sign=field%2?-1.:1.;return 0;
}
