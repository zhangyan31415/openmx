#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "cpu_reclaim.h"
int main(void){
 size_t page=(size_t)sysconf(_SC_PAGESIZE),keep=page+19,total=64*page+31;unsigned char*base=malloc(total+128),*p=base+17;memset(base,0xa7,total+128);
 assert(CpuReclaim_Pages(NULL,0,total)==0);assert(CpuReclaim_Pages(p,total+1,total)==0);assert(CpuReclaim_Pages(p,0,SIZE_MAX)==0);
 size_t n=CpuReclaim_Pages(p,keep,total);if(!n){fprintf(stderr,"FAIL live-prefix page reclaim returned zero\n");return 1;}assert(n%page==0);
 uintptr_t a=((uintptr_t)(p+keep)+page-1)/page*page,b=(uintptr_t)(p+total)/page*page;
 assert(n==b-a);for(size_t i=0;i<keep;i++)assert(p[i]==0xa7);for(size_t i=0;i<17;i++)assert(base[i]==0xa7);for(size_t i=17+total;i<total+128;i++)assert(base[i]==0xa7);
 for(uintptr_t q=a;q<b;q++)assert(*(unsigned char*)q==0);memset(p,0x35,total);for(size_t i=0;i<total;i++)assert(p[i]==0x35);free(base);puts("PASS whole-page reclaim preserves live prefix and allocator-neighbor canaries");return 0;
}
