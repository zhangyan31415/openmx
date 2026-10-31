#ifndef CPU_LOCAL_OLP_H
#define CPU_LOCAL_OLP_H
#include "cpu_local_hs.h"
/* Keep the owner rows and both pre-existing special scratch slots. */
static int CpuLocalOLPPlan(const char*flag,const char*hs,const CpuLocalHSInput*x){
 if(!flag||!strcmp(flag,"0"))return 0;
 if(strcmp(flag,"1"))return -1;
 if(!hs||strcmp(hs,"1")||CpuLocalHSPlan(hs,x)!=1)return -2;
 return 1;
}
#endif
