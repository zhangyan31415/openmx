#ifndef GSOC_INDEX_H
#define GSOC_INDEX_H
#include <stdint.h>
static int gsoc_global(int local,int block,int coord,int procs){return ((local/block)*procs+coord)*block+local%block;}
/* Matches half-open floor(N*r/P) native logical-state intervals. */
static int gsoc_owner(int n,int p,int index){return (int)(((int64_t)(index+1)*p-1)/n);}
#endif
