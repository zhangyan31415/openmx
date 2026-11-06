#ifndef CPU_BAND_OPTIONS_H
#define CPU_BAND_OPTIONS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <errno.h>
#ifndef OPENMX_BAND_PYTHON
#define OPENMX_BAND_PYTHON "python3"
#endif
enum {CPU_BAND_LEGACY=0,CPU_BAND_SLEPC=1,CPU_BAND_GVD=2,CPU_BAND_EVR=3};
typedef struct {int count,solver,band_on,worker;char reference[4096],python[PATH_MAX];} CpuBandConfig;
extern CpuBandConfig CpuBandOptions;
static const char *cpu_band_solver_name(void){
 return CpuBandOptions.solver==CPU_BAND_SLEPC?"slepc":CpuBandOptions.solver==CPU_BAND_GVD?"gvd":CpuBandOptions.solver==CPU_BAND_EVR?"evr":"legacy";
}
#ifdef CPU_BAND_INPUT_IMPLEMENTATION
/* Called while the ordinary native Inputtools file is open. */
static int cpu_band_read_options(void){
 char number[4096],method[4096],*end;long count;int paths=0,explicit_method;
 memset(&CpuBandOptions,0,sizeof(CpuBandOptions));
 if(input_string("Band.NumBands",number,"0")<0)return 1;
 errno=0;count=strtol(number,&end,10);
 if(errno||end==number||*end||count<0||count>INT_MAX)return 1;
 CpuBandOptions.count=(int)count;
 if(input_logical("Band.dispersion",&CpuBandOptions.band_on,0)<0)return 1;
 explicit_method=input_string("Band.Solver",method,"Auto");if(explicit_method<0)return 1;
 if(!strcasecmp(method,"Auto"))CpuBandOptions.solver=count?CPU_BAND_SLEPC:(explicit_method?CPU_BAND_EVR:CPU_BAND_LEGACY);
 else if(!strcasecmp(method,"LAPACK")||!strcasecmp(method,"GVD"))CpuBandOptions.solver=CPU_BAND_GVD;
 else if(!strcasecmp(method,"EVR"))CpuBandOptions.solver=CPU_BAND_EVR;
 else if(!strcasecmp(method,"SLEPc"))CpuBandOptions.solver=CPU_BAND_SLEPC;
 else if(!strcasecmp(method,"Legacy"))CpuBandOptions.solver=CPU_BAND_LEGACY;
 else return 1;
 if(input_string("Band.ReferenceEnergy",CpuBandOptions.reference,"Fermi")<0)return 1;
 if(input_string("Band.SLEPc.Python",CpuBandOptions.python,OPENMX_BAND_PYTHON)<0)return 1;
 if(strcasecmp(CpuBandOptions.reference,"Fermi")){
  errno=0;double energy=strtod(CpuBandOptions.reference,&end);
  if(errno||end==CpuBandOptions.reference||*end||!isfinite(energy))return 1;
 }
 if(count&&!CpuBandOptions.band_on)return 1;
 if(CpuBandOptions.band_on){
  if(input_int("Band.Nkpath",&paths,0)<0||paths<1)return 1;
  if((CpuBandOptions.solver==CPU_BAND_SLEPC&&!count)||(CpuBandOptions.solver==CPU_BAND_LEGACY&&count))return 1;
 }
 CpuBandOptions.worker=CpuBandOptions.band_on&&CpuBandOptions.solver!=CPU_BAND_LEGACY;
 if(CpuBandOptions.worker&&!CpuBandOptions.python[0])return 1;
 return 0;
}
#endif
#endif
