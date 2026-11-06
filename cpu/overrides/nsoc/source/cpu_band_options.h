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
typedef struct {int count;char reference[4096],python[PATH_MAX];} CpuBandConfig;
extern CpuBandConfig CpuBandOptions;
#ifdef CPU_BAND_INPUT_IMPLEMENTATION
/* Called while the ordinary native Inputtools file is open. */
static int cpu_band_read_options(void){
 char number[4096],*end;long count;int dispersion=0,paths=0;
 memset(&CpuBandOptions,0,sizeof(CpuBandOptions));
 if(input_string("Band.NumBands",number,"0")<0)return 1;
 errno=0;count=strtol(number,&end,10);
 if(errno||end==number||*end||count<0||count>INT_MAX)return 1;
 CpuBandOptions.count=(int)count;
 if(input_string("Band.ReferenceEnergy",CpuBandOptions.reference,"Fermi")<0)return 1;
 if(input_string("Band.SLEPc.Python",CpuBandOptions.python,OPENMX_BAND_PYTHON)<0)return 1;
 if(strcasecmp(CpuBandOptions.reference,"Fermi")){
  errno=0;double energy=strtod(CpuBandOptions.reference,&end);
  if(errno||end==CpuBandOptions.reference||*end||!isfinite(energy))return 1;
 }
 if(count){
  if(input_logical("Band.dispersion",&dispersion,0)<0||!dispersion||!CpuBandOptions.python[0])return 1;
  if(input_int("Band.Nkpath",&paths,0)<0||paths<1)return 1;
 }
 return 0;
}
#endif
#endif
