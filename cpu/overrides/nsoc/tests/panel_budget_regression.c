#include <assert.h>
#include <stdio.h>
#include "panel_budget.h"
int main(void){
 /* 1TiB headroom /128 ranks /8 reserve factor =1GiB/rank, exceeding old64MiB. */
 assert(OpenMX_PanelBudget((size_t)1<<40,128)==((size_t)1<<30));
 assert(OpenMX_PanelBudget(0,128)==0);
 assert(OpenMX_PanelBudget(1024,0)==0);
 assert(OpenMX_PanelBudget(1024,2)==64);
 assert(OpenMX_PanelBudget((size_t)-1,1)==(size_t)-1/8);
 puts("panel budget arithmetic PASS");return 0;}
