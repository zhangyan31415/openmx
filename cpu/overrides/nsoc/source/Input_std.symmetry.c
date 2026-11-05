/* CPU input adapter. scf.Symmetry on prepares a checked runtime profile once
 * on rank zero and distributes it before the upstream parser initializes SCF.
 * Previously prepared advanced profiles remain available for developer use. */
#define Input_std CpuOriginal_Input_std
#include "../../../vendor/openmx4.0/source/Input_std.c"
#undef Input_std
#include "cpu_input_bootstrap.h"

static int cpu_input_env_equal(const char *key, const char *value) {
  const char *actual=getenv(key);
  return actual && !strcasecmp(actual,value);
}

void Input_std(char *file) {
  char mode[YOUSO10],group[YOUSO10],policy[YOUSO10],option[YOUSO10],hpolicy[YOUSO10],toggle[YOUSO10];
  double tolerance,covariance;
  int bad=0,any=0,rank=0,enabled=0,has_toggle=0,has_mode=0,hs_fileout=0,soc=0;
  MPI_Comm_rank(mpi_comm_level1,&rank);
  if(!input_open(file))bad=1;
  MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(any){
    if(!rank)fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT: cannot read input on every MPI rank\n");
    MPI_Abort(mpi_comm_level1,95);exit(95);
  }
  if(!bad){
    has_toggle=input_string("scf.Symmetry",toggle,"off");
    enabled=!strcasecmp(toggle,"on");
    if(enabled && (stress_flag || StressConsistency_flag)){
      if(!rank)fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT: stress diagnostics are outside the fixed-geometry symmetry scope\n");
      bad=1;
    }
    if(has_toggle<0 || (strcasecmp(toggle,"on") && strcasecmp(toggle,"off")))bad=1;
    has_mode=input_string("scf.Symmetry.Mode",mode,enabled?"Auto":"Off");
    if(has_mode<0 || (has_toggle>0 && ((!enabled && strcasecmp(mode,"off")) ||
       (enabled && !strcasecmp(mode,"off")))))bad=1;
    input_string("scf.Symmetry.PointGroup",group,"none");
    input_string("scf.Symmetry.HSExport",policy,enabled?"symmetrize":"raw");
    input_string("scf.Symmetry.HamiltonianPolicy",hpolicy,enabled?"Symmetrize":"Free");
    if(input_logical("scf.SpinOrbit.Coupling",&soc,0)<0)bad=1;
    if(input_logical("HS.fileout",&hs_fileout,0)<0)bad=1;
    input_double("scf.Symmetry.Tolerance",&tolerance,enabled?5e-5:1e-5);
    input_double("scf.Symmetry.CovarianceTolerance",&covariance,enabled?1e-5:1e-10);
    if(strcasecmp(mode,"off") && strcasecmp(mode,"auto") && strcasecmp(mode,"user"))bad=1;
    if(strcasecmp(policy,"raw") && strcasecmp(policy,"symmetrize"))bad=1;
    if(strcasecmp(hpolicy,"free") && strcasecmp(hpolicy,"symmetrize"))bad=1;
    if(!strcasecmp(hpolicy,"symmetrize") && (!strcasecmp(mode,"off") || strcasecmp(policy,"symmetrize")))bad=1;
    if(!isfinite(tolerance) || tolerance<=0 || tolerance>.02)bad=1;
    if(!isfinite(covariance) || covariance<1e-14 || covariance>1e-3)bad=1;
    input_string("scf.Symmetry.GeometryPolicy",option,"Strict");if(strcasecmp(option,"Strict"))bad=1;
    input_string("scf.Symmetry.KPointPolicy",option,"Auto");if(strcasecmp(option,"Auto"))bad=1;
    input_string("scf.Symmetry.DensityPolicy",option,"Free");if(strcasecmp(option,"Free"))bad=1;
    input_string("scf.Symmetry.ForcePolicy",option,"Free");if(strcasecmp(option,"Free"))bad=1;
    input_close();
    {
      int selection=enabled+2*soc,minimum,maximum;
      MPI_Allreduce(&selection,&minimum,1,MPI_INT,MPI_MIN,mpi_comm_level1);
      MPI_Allreduce(&selection,&maximum,1,MPI_INT,MPI_MAX,mpi_comm_level1);
      if(minimum!=maximum)bad=1;
    }
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,mpi_comm_level1);
    if(!any && enabled){
      char *environment=NULL,error[PATH_MAX+256]="";
      int length=0,status=0,pinned=cpu_bootstrap_single_cpu_affinity(),all_pinned;
      MPI_Allreduce(&pinned,&all_pinned,1,MPI_INT,MPI_MIN,mpi_comm_level1);
      if(!rank){
        status=cpu_bootstrap_prepare(file,soc?"soc":"nsoc",openmp_threads_num,all_pinned,&environment,error,sizeof(error));
        if(!status)length=(int)strlen(environment)+1;
        else fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT: %s\n",error);
      }
      MPI_Bcast(&status,1,MPI_INT,0,mpi_comm_level1);
      if(status)bad=1;
      else{
        MPI_Bcast(&length,1,MPI_INT,0,mpi_comm_level1);
        if(!rank){if(length<1 || length>(int)CPU_BOOTSTRAP_ENV_MAX+1)bad=1;}
        else{
          if(length<1 || length>(int)CPU_BOOTSTRAP_ENV_MAX+1)bad=1;
          else if(!(environment=malloc((size_t)length)))bad=1;
        }
        MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,mpi_comm_level1);
        if(!any){
          MPI_Bcast(environment,length,MPI_CHAR,0,mpi_comm_level1);
          if(cpu_bootstrap_environment(environment,0,error,sizeof(error)) ||
             cpu_bootstrap_clear_profile(error,sizeof(error)) ||
             cpu_bootstrap_environment(environment,1,error,sizeof(error)))bad=1;
          if(bad)fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT rank=%d: %s\n",rank,error);
        }else bad=1;
      }
      free(environment);
    }else if(any)bad=1;
    if(has_toggle>0 && !enabled && !strcasecmp(mode,"off")){
      char error[256]="";
      if(cpu_bootstrap_clear_profile(error,sizeof(error))){
        fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT rank=%d: %s\n",rank,error);bad=1;
      }
    }
    if(strcasecmp(mode,"off")){
      int legacy=soc && cpu_input_env_equal("OPENMX_CPU_C3","3") && !cpu_input_env_equal("OPENMX_CPU_GSOC","1");
      const char *basis=getenv(legacy?"OPENMX_CPU_C3_BASIS":soc?"OPENMX_CPU_GSOC_BASIS":"OPENMX_CPU_NSOC_BASIS");
      const char *tol=getenv("OPENMX_CPU_INPUT_SYMMETRY_TOL");
      const char *cov=getenv(legacy?"OPENMX_CPU_INPUT_SYMMETRY_COV_TOL":soc?"OPENMX_CPU_GSOC_COV_TOL":"OPENMX_CPU_NSOC_GENERIC_COV_TOL");
      const char *selected=getenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP");
      char *end=NULL;double actual=tol?strtod(tol,&end):NAN;
      if(!cpu_input_env_equal("OPENMX_CPU_INPUT_SYMMETRY_MODE",mode) ||
         !cpu_input_env_equal("OPENMX_CPU_INPUT_SYMMETRY_EXPORT",policy) ||
         !tol || !end || *end || !isfinite(actual) || actual!=tolerance || !selected ||
         (!strcasecmp(mode,"user") && strcasecmp(group,selected)) ||
         (!legacy && !cpu_input_env_equal(soc?"OPENMX_CPU_GSOC":"OPENMX_CPU_NSOC_SYM","1")) ||
         !basis || basis[0]!='/')bad=1;
      end=NULL;actual=cov?strtod(cov,&end):NAN;
      if(!cov||!end||*end||!isfinite(actual)||actual!=covariance)bad=1;
      if(basis){FILE *f=fopen(basis,"rb");if(!f)bad=1;else fclose(f);}
      if(legacy){
        const char *map=getenv("OPENMX_CPU_SOC_NATIVE_MAP");FILE *f=NULL;
        if(!selected || (strcmp(selected,"C3")&&strcmp(selected,"D3")) ||
           strcasecmp(policy,"symmetrize") ||
           !cpu_input_env_equal("OPENMX_CPU_SOC_NATIVE_GROUP",selected?selected:"") ||
           !cpu_input_env_equal("OPENMX_CPU_D3",selected&&!strcmp(selected,"D3")?"1":"0") ||
           !map || map[0]!='/' || !(f=fopen(map,"rb")))bad=1;
        if(f)fclose(f);
        if(selected&&!strcmp(selected,"C3") &&
           (cpu_input_env_equal("OPENMX_CPU_SOC_TR","1")||cpu_input_env_equal("OPENMX_CPU_SOC_REAL","1")))bad=1;
      }
      if(!soc){
        const char *rs=getenv("OPENMX_CPU_REALSPACE_D3");
        int active=rs&&(!strcmp(rs,"scf")||!strcmp(rs,"1"));
        if(active!=!strcasecmp(hpolicy,"symmetrize"))bad=1;
        if(!strcasecmp(policy,"symmetrize")){
          const char *map=getenv("OPENMX_CPU_REALSPACE_MAP");FILE *f=NULL;
          if(!active&&(!rs||strcmp(rs,"export")))bad=1;
          if(!map||map[0]!='/'||!(f=fopen(map,"rb")))bad=1;
          if(f)fclose(f);
        }else if(rs&&strcmp(rs,"off")&&strcmp(rs,"0"))bad=1;
      }
      if(!strcasecmp(policy,"symmetrize") && soc && !legacy &&
         (!cpu_input_env_equal("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE",hs_fileout?"1":"0") ||
          !cpu_input_env_equal("OPENMX_CPU_GSOC_PROJECT","1")))bad=1;
      if(!strcasecmp(policy,"raw") && soc && !legacy &&
         !cpu_input_env_equal("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE","0") &&
         getenv("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE"))bad=1;
    }else if(strcasecmp(policy,"raw"))bad=1;
  }
  MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(any){
    if(!rank)fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT: invalid or unprepared symmetry request; inspect preparation log and input settings\n");
    MPI_Abort(mpi_comm_level1,95);exit(95);
  }
  if(!rank && strcasecmp(mode,"off")){
    printf("CPU_INPUT_SYMMETRY mode=%s selected_group=%s HSExport=%s prepared_profile_consumed=1\n",mode,getenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP"),policy);fflush(stdout);
  }
  CpuOriginal_Input_std(file);
}
