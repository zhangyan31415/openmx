#ifndef CPU_LOCAL_HS_H
#define CPU_LOCAL_HS_H
#include <string.h>
/* This switch removes storage, never atoms, neighbors, integrals or symmetry
 * operations. Geometry changes and unrelated solver paths are excluded. */
typedef struct {
 int md,md_iter,solver,spin,soc,contraction,max_scf,mixing,control_temp,symmetry_mode;
 int hub,constraint,zeeman_s,zeeman_o,core_hole,xmcd;
} CpuLocalHSInput;
static int CpuLocalHSPlan(const char*flag,const CpuLocalHSInput*x){
 if(!flag||!strcmp(flag,"0"))return 0;
 if(strcmp(flag,"1"))return -1;
 if(x->md||x->md_iter!=1||x->solver!=2||x->spin!=3||x->soc!=1||x->contraction||
    x->max_scf!=1||x->mixing!=1||x->control_temp||x->symmetry_mode!=3||
    x->hub||x->constraint||x->zeeman_s||x->zeeman_o||x->core_hole||x->xmcd)return -2;
 return 1;
}
#endif
