/* CPU input adapter. Geometry preparation happens once before MPI launch via
 * workflow.py configure; every native rank verifies that it is consumed.
 * Keep the upstream parser snapshot and unrelated HIP parser untouched. */
#define Input_std CpuOriginal_Input_std
#include "../../../vendor/openmx4.0/source/Input_std.c"
#undef Input_std

static int cpu_input_env_equal(const char *key, const char *value) {
  const char *actual=getenv(key);
  return actual && !strcasecmp(actual,value);
}

void Input_std(char *file) {
  char mode[YOUSO10],group[YOUSO10],policy[YOUSO10],spin[YOUSO10],option[YOUSO10],hpolicy[YOUSO10];
  double tolerance,covariance;
  int bad=0,any=0,rank=0;
  MPI_Comm_rank(mpi_comm_level1,&rank);
  if(!input_open(file))bad=1;
  if(!bad){
    input_string("scf.Symmetry.Mode",mode,"Off");
    input_string("scf.Symmetry.PointGroup",group,"none");
    input_string("scf.Symmetry.HSExport",policy,"raw");
    input_string("scf.Symmetry.HamiltonianPolicy",hpolicy,"Free");
    input_string("scf.SpinOrbit.Coupling",spin,"off");
    input_double("scf.Symmetry.Tolerance",&tolerance,1e-5);
    input_double("scf.Symmetry.CovarianceTolerance",&covariance,1e-10);
    if(strcasecmp(mode,"off") && strcasecmp(mode,"auto") && strcasecmp(mode,"user"))bad=1;
    if(strcasecmp(policy,"raw") && strcasecmp(policy,"symmetrize"))bad=1;
    if(strcasecmp(hpolicy,"free") && strcasecmp(hpolicy,"symmetrize"))bad=1;
    if(!strcasecmp(hpolicy,"symmetrize") && (!strcasecmp(mode,"off") || strcasecmp(policy,"symmetrize")))bad=1;
    if(!isfinite(tolerance) || tolerance<=0 || tolerance>.02)bad=1;
    input_string("scf.Symmetry.GeometryPolicy",option,"Strict");if(strcasecmp(option,"Strict"))bad=1;
    input_string("scf.Symmetry.KPointPolicy",option,"Auto");if(strcasecmp(option,"Auto"))bad=1;
    input_string("scf.Symmetry.DensityPolicy",option,"Free");if(strcasecmp(option,"Free"))bad=1;
    input_string("scf.Symmetry.ForcePolicy",option,"Free");if(strcasecmp(option,"Free"))bad=1;
    if(strcasecmp(mode,"off")){
      const char *basis=getenv(!strcasecmp(spin,"on")?"OPENMX_CPU_GSOC_BASIS":"OPENMX_CPU_NSOC_BASIS");
      const char *tol=getenv("OPENMX_CPU_INPUT_SYMMETRY_TOL");
      const char *cov=getenv(!strcasecmp(spin,"on")?"OPENMX_CPU_GSOC_COV_TOL":"OPENMX_CPU_NSOC_GENERIC_COV_TOL");
      const char *selected=getenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP");
      char *end=NULL;double actual=tol?strtod(tol,&end):NAN;
      if(!cpu_input_env_equal("OPENMX_CPU_INPUT_SYMMETRY_MODE",mode) ||
         !cpu_input_env_equal("OPENMX_CPU_INPUT_SYMMETRY_EXPORT",policy) ||
         !tol || !end || *end || !isfinite(actual) || actual!=tolerance || !selected ||
         (!strcasecmp(mode,"user") && strcasecmp(group,selected)) ||
         !cpu_input_env_equal(!strcasecmp(spin,"on")?"OPENMX_CPU_GSOC":"OPENMX_CPU_NSOC_SYM","1") ||
         !basis || basis[0]!='/')bad=1;
      end=NULL;actual=cov?strtod(cov,&end):NAN;
      if(!cov||!end||*end||!isfinite(actual)||actual!=covariance)bad=1;
      if(basis){FILE *f=fopen(basis,"rb");if(!f)bad=1;else fclose(f);}
      if(strcasecmp(spin,"on")){
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
      if(!strcasecmp(policy,"symmetrize") && !strcasecmp(spin,"on") &&
         (!cpu_input_env_equal("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE","1") ||
          !cpu_input_env_equal("OPENMX_CPU_GSOC_PROJECT","1")))bad=1;
      if(!strcasecmp(policy,"raw") && !strcasecmp(spin,"on") &&
         !cpu_input_env_equal("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE","0") &&
         getenv("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE"))bad=1;
    }else if(strcasecmp(policy,"raw"))bad=1;
    input_close();
  }
  MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(any){
    if(!rank)fprintf(stderr,"CPU_INPUT_SYMMETRY_REJECT: unsupported or unprepared input symmetry request; use cpu/tools/workflow.py configure, then run its complete configuration\n");
    MPI_Abort(mpi_comm_level1,95);exit(95);
  }
  if(!rank && strcasecmp(mode,"off")){
    printf("CPU_INPUT_SYMMETRY mode=%s selected_group=%s HSExport=%s prepared_profile_consumed=1\n",mode,getenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP"),policy);fflush(stdout);
  }
  CpuOriginal_Input_std(file);
}
