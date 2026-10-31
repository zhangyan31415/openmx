/* Shared by the solver and its isolated numerical configuration test. */
#include <math.h>
#include <stdlib.h>
int CpuC3_GeometryToleranceAng(const char *text,double *value)
{
  char *end;
  double parsed;
  if(!text||!text[0]||!value)return 0;
  parsed=strtod(text,&end);
  if(*end||!isfinite(parsed)||parsed<=0.0||parsed>0.02)return 0;
  *value=parsed;
  return 1;
}
double CpuC3_PeriodicDistanceBohr(const double *dr,const double *cell)
{
  double cross[3][3],det,frac[3],wrapped[3]={0.,0.,0.},norm=0.;
  int i,j;
  for(i=0;i<3;i++){
    const double *a=cell+3*((i+1)%3),*b=cell+3*((i+2)%3);
    cross[i][0]=a[1]*b[2]-a[2]*b[1];
    cross[i][1]=a[2]*b[0]-a[0]*b[2];
    cross[i][2]=a[0]*b[1]-a[1]*b[0];
  }
  det=cell[0]*cross[0][0]+cell[1]*cross[0][1]+cell[2]*cross[0][2];
  if(!isfinite(det)||fabs(det)<1e-12)return NAN;
  for(i=0;i<3;i++){
    frac[i]=0.;for(j=0;j<3;j++)frac[i]+=cross[i][j]*dr[j]/det;
    if(!isfinite(frac[i]))return NAN;
    frac[i]-=nearbyint(frac[i]);
  }
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)wrapped[i]+=frac[j]*cell[3*j+i];
    norm+=wrapped[i]*wrapped[i];
  }
  return sqrt(norm);
}
