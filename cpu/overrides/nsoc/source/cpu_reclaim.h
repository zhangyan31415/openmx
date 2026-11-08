#ifndef OPENMX_CPU_RECLAIM_H
#define OPENMX_CPU_RECLAIM_H
#include <stddef.h>
#include <stdint.h>
#if defined(__linux__)
#include <sys/mman.h>
#include <unistd.h>
#endif
/* Return only complete pages strictly inside dead allocation payload. Retain
 * the pointer/virtual capacity: a future explicit rebuild can write it again.
 * Boundary pages can contain live prefix or allocator/neighbor data. */
static size_t CpuReclaim_Pages(void *allocation,size_t live_bytes,size_t total_bytes){
#if defined(__linux__) && defined(MADV_DONTNEED)
 uintptr_t begin=(uintptr_t)allocation,first,last;long raw_page=sysconf(_SC_PAGESIZE);
 if(!allocation||live_bytes>total_bytes||raw_page<=0||total_bytes>UINTPTR_MAX-begin)return 0;
 size_t page=(size_t)raw_page;first=begin+live_bytes;last=begin+total_bytes;
 size_t rem=first%page;if(rem){size_t add=page-rem;if(first>UINTPTR_MAX-add)return 0;first+=add;}
 last-=last%page;if(first>=last)return 0;
 size_t bytes=(size_t)(last-first);return madvise((void*)first,bytes,MADV_DONTNEED)==0?bytes:0;
#else
 (void)allocation;(void)live_bytes;(void)total_bytes;return 0;
#endif
}
#endif
