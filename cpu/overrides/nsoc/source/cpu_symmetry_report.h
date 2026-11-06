/* Presentation only: preserve the rank-zero preparation report in native .out. */
#ifndef OPENMX_CPU_SYMMETRY_REPORT_H
#define OPENMX_CPU_SYMMETRY_REPORT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static int cpu_symmetry_report_copy(FILE *destination) {
  const char *path=getenv("OPENMX_CPU_INPUT_SYMMETRY_REPORT");
  char buffer[4096];size_t count;int bad=0;FILE *source;
  if(!path || !*path)return 0;
  source=fopen(path,"rb");
  if(!source)return 1;
  while((count=fread(buffer,1,sizeof(buffer),source)))
    if(fwrite(buffer,1,count,destination)!=count){bad=1;break;}
  if(ferror(source))bad=1;
  if(fclose(source))bad=1;
  return bad;
}

static int cpu_symmetry_report_initial(const char *directory,const char *name) {
  char path[PATH_MAX];FILE *out;int bad;
  if(!getenv("OPENMX_CPU_INPUT_SYMMETRY_REPORT"))return 0;
  if(snprintf(path,sizeof(path),"%s%s.out",directory,name)>=(int)sizeof(path))return 1;
  out=fopen(path,"w");if(!out)return 1;
  fputs("OpenMX: input initialized; calculation in progress.\n"
        "Final input and physical results are merged here at normal completion.\n\n",out);
  bad=cpu_symmetry_report_copy(out);
  if(fclose(out))bad=1;
  return bad;
}
#endif
