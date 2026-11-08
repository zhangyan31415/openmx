#ifndef OPENMX_PANEL_BUDGET_H
#define OPENMX_PANEL_BUDGET_H
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
/* At most one eighth of current free memory is available to the new planner,
 * shared equally by this node's ranks. Divide first to avoid multiplication overflow. */
static size_t OpenMX_PanelBudget(size_t free_bytes,int local_ranks){
 return local_ranks>0?(free_bytes/8)/(size_t)local_ranks:0;
}
static size_t OpenMX_PanelHeadroom(void){
 FILE*f=fopen("/proc/meminfo","r");char line[4096];unsigned long long kb=0;size_t headroom=0;
 if(f){while(fgets(line,sizeof(line),f))if(sscanf(line,"MemAvailable: %llu kB",&kb)==1){headroom=kb>SIZE_MAX/1024?SIZE_MAX:(size_t)kb*1024;break;}fclose(f);}
 /* A node can have more RAM than this task's cgroup allocation. Respect every
  * numeric cgroup-v2 ancestor limit; 'max' is not a numeric limit. */
 f=fopen("/proc/self/cgroup","r");char cg[4096]={0};
 if(f){while(fgets(line,sizeof(line),f))if(!strncmp(line,"0::/",4)){size_t n=strcspn(line+3,"\r\n");if(n+strlen("/sys/fs/cgroup")+1<sizeof(cg)){strcpy(cg,"/sys/fs/cgroup");strncat(cg,line+3,n);}break;}fclose(f);}
 while(cg[0]){
  char path[4200];unsigned long long limit=0,current=0;int has_limit=0,has_current=0;
  snprintf(path,sizeof(path),"%s/memory.max",cg);f=fopen(path,"r");if(f){has_limit=fscanf(f,"%llu",&limit)==1;fclose(f);}
  snprintf(path,sizeof(path),"%s/memory.current",cg);f=fopen(path,"r");if(f){has_current=fscanf(f,"%llu",&current)==1;fclose(f);}
  if(has_limit&&has_current){unsigned long long free_bytes=current<limit?limit-current:0;size_t available=free_bytes>SIZE_MAX?SIZE_MAX:(size_t)free_bytes;if(available<headroom)headroom=available;}
  if(!strcmp(cg,"/sys/fs/cgroup")||strlen(cg)<=strlen("/sys/fs/cgroup"))break;
  char*slash=strrchr(cg,'/');if(!slash)break;*slash=0;
 }
 return headroom;
}
#endif
