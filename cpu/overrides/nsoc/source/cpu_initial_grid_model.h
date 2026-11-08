#ifndef CPU_INITIAL_GRID_MODEL_H
#define CPU_INITIAL_GRID_MODEL_H
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int n[3],frame_count;
  double origin[3],inverse[9],*frames,*values;
  int warp_count,*wave;
  double warp_inverse[4],*warp_coeff;
  int isolated_z;
} CpuInitialGridModel;

static void CpuInitialGridFree(CpuInitialGridModel *m)
{free(m->frames);free(m->values);free(m->wave);free(m->warp_coeff);memset(m,0,sizeof(*m));}

static int CpuInitialGridRead(const char *path,CpuInitialGridModel *m)
{
  FILE *f=fopen(path,"rb");char magic[8];int dims[4];size_t count=1;
  memset(m,0,sizeof(*m));if(!f)return -1;
  if(sizeof(int)!=4||fread(magic,1,8,f)!=8||
     (memcmp(magic,"OMXIGR1",8)&&memcmp(magic,"OMXIGR2",8))||
     fread(dims,sizeof(int),4,f)!=4)goto bad;
  for(int i=0;i<3;i++){
    if(dims[i]<1||count>SIZE_MAX/(size_t)dims[i])goto bad;
    m->n[i]=dims[i];count*=dims[i];
  }
  if(dims[3]<1||(size_t)dims[3]>SIZE_MAX/(4*sizeof(double))||count>SIZE_MAX/sizeof(double))goto bad;
  m->frame_count=dims[3];
  m->frames=malloc((size_t)m->frame_count*4*sizeof(double));
  m->values=malloc(count*sizeof(double));
  if(!m->frames||!m->values||fread(m->origin,sizeof(double),3,f)!=3||
     fread(m->inverse,sizeof(double),9,f)!=9||
     fread(m->frames,sizeof(double),(size_t)m->frame_count*4,f)!=(size_t)m->frame_count*4)goto bad;
  if(!memcmp(magic,"OMXIGR2",8)){
    m->isolated_z=1;
    if(fread(&m->warp_count,sizeof(int),1,f)!=1||m->warp_count<1||m->warp_count>INT32_MAX/2||
       (size_t)m->frame_count>SIZE_MAX/((size_t)m->warp_count*6*sizeof(double)))goto bad;
    size_t nc=(size_t)m->frame_count*m->warp_count*6;
    m->wave=malloc((size_t)m->warp_count*2*sizeof(int));m->warp_coeff=malloc(nc*sizeof(double));
    if(!m->wave||!m->warp_coeff||fread(m->warp_inverse,sizeof(double),4,f)!=4||
       fread(m->wave,sizeof(int),(size_t)m->warp_count*2,f)!=(size_t)m->warp_count*2||
       fread(m->warp_coeff,sizeof(double),nc,f)!=nc)goto bad;
    for(int i=0;i<4;i++)if(!isfinite(m->warp_inverse[i]))goto bad;
    for(size_t i=0;i<nc;i++)if(!isfinite(m->warp_coeff[i]))goto bad;
  }
  if(fread(m->values,sizeof(double),count,f)!=count||fgetc(f)!=EOF)goto bad;
  for(int i=0;i<3;i++)if(!isfinite(m->origin[i]))goto bad;
  for(int i=0;i<9;i++)if(!isfinite(m->inverse[i]))goto bad;
  for(size_t i=0;i<(size_t)m->frame_count*4;i++)if(!isfinite(m->frames[i]))goto bad;
  for(size_t i=0;i<count;i++)if(!isfinite(m->values[i]))goto bad;
  fclose(f);return 0;
bad:
  fclose(f);CpuInitialGridFree(m);return -1;
}

/* Periodic trilinear interpolation is only an initial-guess operation.
 * It never replaces the native quadrature of a final SCF iteration. */
static double CpuInitialGridSample(const CpuInitialGridModel *m,const double r[3])
{
  double value=0.;
  for(int frame=0;frame<m->frame_count;frame++){
    const double *t=m->frames+4*(size_t)frame;
    double u[3]={0.,0.,0.};
    if(m->warp_count){
      double f0=m->warp_inverse[0]*r[0]+m->warp_inverse[1]*r[1];
      double f1=m->warp_inverse[2]*r[0]+m->warp_inverse[3]*r[1];
      for(int wv=0;wv<m->warp_count;wv++){
        double phase=6.2831853071795864769*(m->wave[2*wv]*f0+m->wave[2*wv+1]*f1);
        double co=cos(phase),si=sin(phase);
        const double *v=m->warp_coeff+((size_t)frame*m->warp_count+wv)*6;
        for(int a=0;a<3;a++)u[a]+=v[2*a]*co+v[2*a+1]*si;
      }
    }
    double c=cos(t[0]),s=sin(t[0]),dx=r[0]-u[0]-t[1],dy=r[1]-u[1]-t[2];
    double q[3]={c*dx+s*dy-m->origin[0],-s*dx+c*dy-m->origin[1],r[2]-u[2]-t[3]-m->origin[2]};
    int i0[3],outside=0;double w[3];
    for(int a=0;a<3;a++){
      double f=0.;for(int b=0;b<3;b++)f+=m->inverse[3*a+b]*q[b];
      if(!isfinite(f))return NAN;
      if(a==2&&m->isolated_z&&(f<0.||f>=1.)){outside=1;break;}
      f-=floor(f);f*=m->n[a];
      i0[a]=(int)floor(f);if(i0[a]>=m->n[a])i0[a]=0;
      w[a]=f-floor(f);
    }
    if(outside)continue;
    for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int k=0;k<2;k++){
      size_t i=(size_t)((i0[0]+a)%m->n[0]);
      size_t j=(size_t)((i0[1]+b)%m->n[1]);
      size_t l=(size_t)((i0[2]+k)%m->n[2]);
      value+=(a?w[0]:1.-w[0])*(b?w[1]:1.-w[1])*(k?w[2]:1.-w[2])*
        m->values[(i*(size_t)m->n[1]+j)*(size_t)m->n[2]+l];
    }
  }
  return value;
}
#endif
