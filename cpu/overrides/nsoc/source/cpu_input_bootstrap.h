/* Rank-zero native input preparation. No shell evaluation or run deadline. */
#ifndef OPENMX_CPU_INPUT_BOOTSTRAP_H
#define OPENMX_CPU_INPUT_BOOTSTRAP_H
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#ifndef OPENMX_SYMMETRY_BOOTSTRAP
#define OPENMX_SYMMETRY_BOOTSTRAP ""
#endif
#ifndef OPENMX_SYMMETRY_PYTHON
#define OPENMX_SYMMETRY_PYTHON "python3"
#endif
#define CPU_BOOTSTRAP_ENV_MAX (16U*1024U*1024U)
extern char **environ;

/* Automatic native preparation supplies a complete profile. A stale developer
 * selector inherited from an earlier run must not silently add another path. */
static int cpu_bootstrap_clear_profile(char *error,size_t error_size) {
  int i=0;
  while(environ[i]){
    if(!strncmp(environ[i],"OPENMX_CPU_",11) &&
       strncmp(environ[i],"OPENMX_CPU_SYMMETRY_HELPER=",27) &&
       strncmp(environ[i],"OPENMX_CPU_SYMMETRY_PYTHON=",27)){
      char *key=strdup(environ[i]),*end;
      if(!key){snprintf(error,error_size,"profile cleanup allocation failed");return 1;}
      end=strchr(key,'=');if(end)*end='\0';
      if(unsetenv(key)){free(key);snprintf(error,error_size,"cannot clear inherited CPU selector");return 1;}
      free(key);
    }else i++;
  }
  return 0;
}

/* Validate the entire protocol before changing this rank's environment. The
 * prefix is a capability boundary: helpers cannot change MPI/loader settings. */
static int cpu_bootstrap_environment(const char *text, int apply,
                                     char *error, size_t error_size) {
  char *copy=strdup(text), *p, *next, *eq, *names[4096];
  int count=0, i, result=1;
  if(!copy){snprintf(error,error_size,"environment allocation failed");return 1;}
  for(p=copy;*p;p=next){
    next=strchr(p,'\n');
    if(next)*next++='\0';else next=p+strlen(p);
    eq=strchr(p,'=');
    if(!eq || eq==p || count==4096)goto invalid;
    *eq='\0';
    if(strncmp(p,"OPENMX_CPU_",11) || strlen(p)>255)goto invalid;
    for(i=0;p[i];i++)if(!((p[i]>='A'&&p[i]<='Z') ||
         (p[i]>='a'&&p[i]<='z') || (p[i]>='0'&&p[i]<='9') || p[i]=='_'))goto invalid;
    for(i=0;i<count;i++)if(!strcmp(names[i],p))goto invalid;
    names[count++]=p;
    for(i=1;eq[i];i++)if((unsigned char)eq[i]<32 || (unsigned char)eq[i]==127)goto invalid;
  }
  if(!count)goto invalid;
  if(apply)for(i=0;i<count;i++){
    p=names[i];eq=p+strlen(p)+1;
    if(setenv(p,eq,1)){snprintf(error,error_size,"cannot apply helper environment: %s",strerror(errno));goto done;}
  }
  result=0;goto done;
invalid:
  snprintf(error,error_size,"invalid, duplicate, or unsupported helper environment entry");
done:
  free(copy);return result;
}

static const char *cpu_bootstrap_helper(char *resolved,size_t capacity) {
  const char *override=getenv("OPENMX_CPU_SYMMETRY_HELPER");
  char executable[PATH_MAX], *slash;
  ssize_t length;
  if(override && *override)return override;
  length=readlink("/proc/self/exe",executable,sizeof(executable)-1);
  if(length>0){
    executable[length]='\0';slash=strrchr(executable,'/');
    if(slash){
      *slash='\0';
      if(snprintf(resolved,capacity,"%s/resources/cpu/geometry/bootstrap.py",executable)<(int)capacity && access(resolved,R_OK)==0)return resolved;
      if(snprintf(resolved,capacity,"%s/../share/openmx/cpu/geometry/bootstrap.py",executable)<(int)capacity && access(resolved,R_OK)==0)return resolved;
    }
  }
  return OPENMX_SYMMETRY_BOOTSTRAP;
}

static int cpu_bootstrap_single_cpu_affinity(void) {
  char line[4096],extra;int cpu,one=0;FILE *f=fopen("/proc/self/status","r");
  if(!f)return 0;
  while(fgets(line,sizeof(line),f))if(!strncmp(line,"Cpus_allowed_list:",18)){
    one=sscanf(line+18,"%d %c",&cpu,&extra)==1 && cpu>=0;break;
  }
  fclose(f);return one;
}

static int cpu_bootstrap_prepare(const char *input,const char *stage,int threads,int ranks_pinned,
                                 char **environment,char *error,size_t error_size) {
  char absolute[PATH_MAX], cwd[PATH_MAX], directory[PATH_MAX];
  char protocol[PATH_MAX+32], logfile[PATH_MAX+32], helper_path[PATH_MAX];
  const char *helper=cpu_bootstrap_helper(helper_path,sizeof(helper_path));
  const char *python=getenv("OPENMX_CPU_SYMMETRY_PYTHON");
  char *argv[16],thread_value[32];
  posix_spawn_file_actions_t actions;
  pid_t child;
  int rc,status;
  FILE *f=NULL;
  struct stat st;
  size_t length;
  *environment=NULL;
  if(!python || !*python)python=OPENMX_SYMMETRY_PYTHON;
  if(!helper || !*helper){snprintf(error,error_size,"symmetry preparation helper is not installed");return 1;}
  if(!realpath(input,absolute) || !getcwd(cwd,sizeof(cwd))){
    snprintf(error,error_size,"cannot resolve input or launch directory: %s",strerror(errno));return 1;
  }
  if(snprintf(directory,sizeof(directory),"%s/.openmx-symmetry-XXXXXX",cwd)>=(int)sizeof(directory) ||
     strlen(directory)+32>=sizeof(protocol)){
    snprintf(error,error_size,"launch directory path is too long");return 1;
  }
  if(!mkdtemp(directory)){snprintf(error,error_size,"cannot create symmetry preparation directory: %s",strerror(errno));return 1;}
  snprintf(protocol,sizeof(protocol),"%s/environment.txt",directory);
  snprintf(logfile,sizeof(logfile),"%s/bootstrap.log",directory);
  printf("CPU_INPUT_SYMMETRY_PREPARE stage=%s directory=%s\n",stage,directory);fflush(stdout);
  argv[0]=(char*)python;argv[1]=(char*)helper;argv[2]="--input";argv[3]=absolute;
  argv[4]="--output";argv[5]=directory;argv[6]="--stage";argv[7]=(char*)stage;
  snprintf(thread_value,sizeof(thread_value),"%d",threads);
  argv[8]="--environment-file";argv[9]=protocol;
  argv[10]="--threads";argv[11]=thread_value;
  argv[12]="--rank-pinned";argv[13]=ranks_pinned?"1":"0";argv[14]=NULL;
  rc=posix_spawn_file_actions_init(&actions);
  if(rc){snprintf(error,error_size,"cannot initialize preparation log: %s",strerror(rc));return 1;}
  if(!rc)rc=posix_spawn_file_actions_addopen(&actions,STDOUT_FILENO,logfile,O_WRONLY|O_CREAT|O_EXCL,0600);
  if(!rc)rc=posix_spawn_file_actions_adddup2(&actions,STDOUT_FILENO,STDERR_FILENO);
  if(rc){posix_spawn_file_actions_destroy(&actions);snprintf(error,error_size,"cannot configure preparation log: %s",strerror(rc));return 1;}
  rc=posix_spawnp(&child,python,&actions,NULL,argv,environ);
  posix_spawn_file_actions_destroy(&actions);
  if(rc){snprintf(error,error_size,"cannot start symmetry preparation: %s (log %s)",strerror(rc),logfile);return 1;}
  do{rc=waitpid(child,&status,0);}while(rc<0 && errno==EINTR);
  if(rc<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=0){
    snprintf(error,error_size,"symmetry preparation failed; inspect %s",logfile);return 1;
  }
  f=fopen(protocol,"rb");
  if(!f || fstat(fileno(f),&st) || !S_ISREG(st.st_mode) || st.st_size<=0 || st.st_size>CPU_BOOTSTRAP_ENV_MAX){
    if(f)fclose(f);
    snprintf(error,error_size,"missing or invalid preparation environment; inspect %s",logfile);return 1;
  }
  length=(size_t)st.st_size;
  *environment=malloc(length+1);
  if(!*environment){fclose(f);snprintf(error,error_size,"environment allocation failed");return 1;}
  if(fread(*environment,1,length,f)!=length || ferror(f) || memchr(*environment,'\0',length)){
    fclose(f);free(*environment);*environment=NULL;
    snprintf(error,error_size,"cannot read preparation environment");return 1;
  }
  fclose(f);(*environment)[length]='\0';
  if(cpu_bootstrap_environment(*environment,0,error,error_size)){
    free(*environment);*environment=NULL;return 1;
  }
  /* Detailed machine diagnostics remain in bootstrap.log. The input adapter
   * prints the readable report after the native filename/path are known. */
  return 0;
}
#endif
