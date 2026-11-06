/* Isolate PETSc's MPI/runtime from the parent OpenMX MPI ABI. No shell or timer. */
#ifndef CPU_BAND_BRIDGE_H
#define CPU_BAND_BRIDGE_H
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#ifndef OPENMX_BAND_HELPER
#define OPENMX_BAND_HELPER ""
#endif
extern char **environ;
static const char *cpu_band_helper(char *path,size_t size){
 char executable[PATH_MAX];ssize_t n=readlink("/proc/self/exe",executable,sizeof(executable)-1);
 if(n>0&&n<(ssize_t)sizeof(executable)-1){
  executable[n]='\0';char *slash=strrchr(executable,'/');if(slash)*slash='\0';
  if(snprintf(path,size,"%s/resources/cpu/tools/band_solver.py",executable)<(int)size&&access(path,R_OK)==0)return path;
  if(snprintf(path,size,"%s/../share/openmx/cpu/tools/band_solver.py",executable)<(int)size&&access(path,R_OK)==0)return path;
 }
 return OPENMX_BAND_HELPER;
}
static int cpu_band_environment_keep(const char *entry){
 const char *prefix[]={"PMI_","PMIX_","I_MPI_","OMPI_","MPI_","HYDRA_","UCX_","FI_","KMP_","LD_PRELOAD=","LD_LIBRARY_PATH=",NULL};
 for(int i=0;prefix[i];i++)if(!strncmp(entry,prefix[i],strlen(prefix[i])))return 0;
 return 1;
}
static int cpu_band_spawn(char *const *argv,const char *log){
 size_t n=0;while(environ[n])n++;
 char **environment=calloc(n+1,sizeof(char*));if(!environment)return 1;
 size_t j=0;for(size_t i=0;i<n;i++)if(cpu_band_environment_keep(environ[i]))environment[j++]=environ[i];
 posix_spawn_file_actions_t actions;int rc=posix_spawn_file_actions_init(&actions),status=0;pid_t child;
 if(rc){free(environment);return 1;}
 if(log){
  rc=posix_spawn_file_actions_addopen(&actions,STDOUT_FILENO,log,O_WRONLY|O_CREAT|O_EXCL,0600);
  if(!rc)rc=posix_spawn_file_actions_adddup2(&actions,STDOUT_FILENO,STDERR_FILENO);
 }
 if(!rc)rc=posix_spawnp(&child,argv[0],&actions,NULL,argv,environment);
 posix_spawn_file_actions_destroy(&actions);free(environment);
 if(rc)return 1;
 do{rc=waitpid(child,&status,0);}while(rc<0&&errno==EINTR);
 return rc<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0;
}
static int cpu_band_check_backend(const char *python){
 char path[PATH_MAX];const char *helper=cpu_band_helper(path,sizeof(path));
 if(!helper||!*helper)return 1;
 char *args[]={(char*)python,(char*)helper,"--check",NULL};
 return cpu_band_spawn(args,NULL);
}
#endif
