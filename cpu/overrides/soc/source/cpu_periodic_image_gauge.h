#ifndef CPU_PERIODIC_IMAGE_GAUGE_H
#define CPU_PERIODIC_IMAGE_GAUGE_H
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/* Metadata may use different lattice representatives than native OpenMX.
 * If f_ref(a)=f_native(a)+m(a), the native image is
 * s_native(g,a)=s_ref(g,a)+m(p_g(a))-W_g m(a).
 * Only exact lattice translations within the fixed-geometry tolerance are
 * accepted. No atomic positions or physical operators are projected here. */
static int cpu_periodic_image_rebase(int nat,int ng,const double cell[9],
 double *reference,double **native,double scale,int **perm,int w[][9],
 int **image,double tolerance,int *moved,char *error,size_t error_size)
{
 int a,i,j,g,changed=0; int *offset=NULL;
 double cross[9],det;
#define CPU_GAUGE_FAIL(message) do {snprintf(error,error_size,"%s",message);free(offset);return -1;} while(0)
 *moved=0;
 if(nat<1||ng<1||!isfinite(scale)||!isfinite(tolerance)||tolerance<=0)
   CPU_GAUGE_FAIL("invalid periodic image rebase dimensions or units");
 if((size_t)nat>SIZE_MAX/(3*sizeof(int)))CPU_GAUGE_FAIL("periodic image scratch overflow");
 for(i=0;i<9;i++)if(!isfinite(cell[i]))CPU_GAUGE_FAIL("nonfinite periodic image cell");
 for(i=0;i<3;i++){
   int next=(i+1)%3,last=(i+2)%3;
   for(j=0;j<3;j++)cross[3*i+j]=cell[3*next+(j+1)%3]*cell[3*last+(j+2)%3]-cell[3*next+(j+2)%3]*cell[3*last+(j+1)%3];
 }
 det=cell[0]*cross[0]+cell[1]*cross[1]+cell[2]*cross[2];
 if(!isfinite(det)||fabs(det)<1e-20)CPU_GAUGE_FAIL("singular periodic image cell");
 offset=calloc((size_t)nat*3,sizeof(int));
 if(!offset)CPU_GAUGE_FAIL("periodic image scratch allocation failed");
 for(a=0;a<nat;a++){
   double delta[3];int shifted=0;
   for(i=0;i<3;i++){
     double x=native[a+1][i+1]*scale;
     if(!isfinite(x)||!isfinite(reference[3*a+i]))CPU_GAUGE_FAIL("nonfinite periodic image coordinates");
     delta[i]=reference[3*a+i]-x;
   }
   for(i=0;i<3;i++){
     double f=0,integer;
     for(j=0;j<3;j++)f+=delta[j]*cross[3*i+j]/det;
     integer=round(f);
     if(!isfinite(integer)||integer<INT_MIN||integer>INT_MAX)CPU_GAUGE_FAIL("periodic image offset overflow");
     offset[3*a+i]=(int)integer;shifted|=(integer!=0);
   }
   for(i=0;i<3;i++){
     double residual=delta[i];
     for(j=0;j<3;j++)residual-=offset[3*a+j]*cell[3*j+i];
     if(fabs(residual)>tolerance)CPU_GAUGE_FAIL("geometry differs from metadata: not a lattice translation");
   }
   changed+=shifted;
 }
 /* Validate every index and widened sum before mutating any metadata. */
 for(g=0;g<ng;g++)for(a=0;a<nat;a++){
   int p=perm[g][a];
   if(p<0||p>=nat)CPU_GAUGE_FAIL("invalid periodic image permutation");
   for(i=0;i<3;i++){
     long double value=(long double)image[g][3*a+i]+offset[3*p+i];
     for(j=0;j<3;j++)value-=(long double)w[g][3*i+j]*offset[3*a+j];
     if(value<INT_MIN||value>INT_MAX)CPU_GAUGE_FAIL("rebased periodic image overflow");
   }
 }
 for(g=0;g<ng;g++)for(a=0;a<nat;a++)for(i=0;i<3;i++){
   long double value=(long double)image[g][3*a+i]+offset[3*perm[g][a]+i];
   for(j=0;j<3;j++)value-=(long double)w[g][3*i+j]*offset[3*a+j];
   image[g][3*a+i]=(int)value;
 }
 for(a=0;a<nat;a++)for(i=0;i<3;i++)reference[3*a+i]=native[a+1][i+1]*scale;
 free(offset);*moved=changed;return 0;
#undef CPU_GAUGE_FAIL
}
#endif
