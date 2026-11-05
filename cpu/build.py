#!/usr/bin/env python3
"""Compile source-only CPU OpenMX + bundled ELPA into a new output directory."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,os,shlex,socket,subprocess,time,shutil,sys
REPO=Path(__file__).resolve().parents[1]
BASE=REPO/'cpu/overrides/soc/SOURCES.json'
DEFAULT_ONEAPI=os.environ.get('ONEAPI_ROOT','/opt/intel/oneapi')

def runtime_resource_files():
 files=[]
 for directory in ('cpu/geometry','cpu/overrides/nsoc/metadata'):
  files.extend(p for p in (REPO/directory).rglob('*')if p.is_file()and p.suffix in ('.py','.json'))
 files.extend(REPO/'cpu/tools'/name for name in ('workflow.py','workflow_node.py','input_symmetry.py'))
 return sorted(set(files))

def copy_runtime_resources(output):
 root=Path(output)/'resources'
 for source in runtime_resource_files():
  dest=root/source.relative_to(REPO);dest.parent.mkdir(parents=True,exist_ok=True)
  shutil.copy2(source,dest)
 return root
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def local(path):
 p=Path(path)
 if p.is_absolute()or'..'in p.parts:raise ValueError('source/include must be a repository-relative path: '+str(path))
 resolved=(REPO/p).resolve()
 if not resolved.is_relative_to(REPO):raise ValueError('source escapes repository: '+str(path))
 if not resolved.exists():raise ValueError('missing source/include: '+str(path))
 return resolved

def source_plan(selection=None, profile="soc"):
 plan=json.loads(BASE.read_text());sel=selection or{}
 if plan['unresolved_objects']or plan['unresolved_includes']:raise ValueError('source closure incomplete')
 units=plan['units']
 if profile=='stock':
  stock=[]
  for u in units:
   if '/vendor/'in u['source']:stock.append(u);continue
   name={'Cluster.compact.o':'Cluster_DFT_NonCol'}.get(u['object'],u['object'].split('.')[0])
   native='cpu/vendor/openmx4.0/source/'+name+'.c'
   if (REPO/native).is_file():stock.append({'object':Path(native).stem+'.o','source':native,'language':'c','flags':[],'historical_object':u['historical_object']})
  units=stock;plan['units']=units;plan['include_dirs']=['cpu/vendor/openmx4.0/source','cpu/vendor/openmx4.0/source/elpa-2018.05.001']
 known={u['object']for u in units}
 for obj,src in sel.get('replacements',{}).items():
  if obj not in known:raise ValueError('unknown replacement object: '+obj)
  local(src)
  for u in units:
   if u['object']==obj:u['source']=src
 units.extend(sel.get('extra_units',[]));seen=set()
 for u in units:
  if Path(u['object']).name!=u['object']or not u['object'].endswith('.o')or u['object']in seen:raise ValueError('invalid/duplicate object: '+u['object'])
  seen.add(u['object']);local(u['source'])
  if u['language']not in ['c','fortran']:raise ValueError('unknown language')
  for flag in u.get('flags',[]):
   if flag.startswith('-I')and len(flag)>2:local(flag[2:])
 plan['include_dirs']=sel.get('include_dirs',[])+plan['include_dirs']
 for d in plan['include_dirs']:local(d)
 plan['capabilities']={'nsoc_symmetry':any(u['object']=='cpu_nsoc_sym.o'for u in units),'legacy_soc':profile!='stock','generic_soc':profile!='stock','compact_soc':profile!='stock','multinode':True,**sel.get('capabilities',{})}
 return plan

def command_plan(plan,out,oneapi,cpu_target):
 mpi=oneapi/'mpi/2021.7.1/bin';mkl=oneapi/'mkl/2022.2.1'
 common=['-O3','-x'+cpu_target,'-ip','-no-prec-div','-qopenmp']
 incs=['-I'+str(local(p))for p in plan['include_dirs']]+['-I'+str(mkl/'include'),'-I'+str(mkl/'include/fftw')]
 commands=[]
 for u in plan['units']:
  flags=[];raw=u.get('flags',[]);i=0
  while i<len(raw):
   f=raw[i]
   if f=='-include':i+=1;flags += ['-include',str(local(raw[i]))]
   elif f.startswith('-I'):flags += ['-I'+str(local(f[2:]))]
   else:flags.append(f)
   i+=1
  cmd=[str(mpi/('mpiifort'if u['language']=='fortran'else'mpiicc')),*common,*incs,*flags,'-c',str(local(u['source'])),'-o',str(out/'objects'/u['object'])]
  if u['source']=='cpu/overrides/nsoc/source/Input_std.symmetry.c':
   cmd[1:1]=['-DOPENMX_SYMMETRY_BOOTSTRAP="'+str(out/'resources/cpu/geometry/bootstrap.py')+'"',
             '-DOPENMX_SYMMETRY_PYTHON="'+sys.executable+'"']
  commands.append({'label':u['object'],'language':u['language'],'command':cmd})
 libs=['-L'+str(mkl/'lib/intel64'),'-lmkl_scalapack_lp64','-lmkl_intel_lp64','-lmkl_intel_thread','-lmkl_core','-lifcore','-lmkl_blacs_intelmpi_lp64','-liomp5','-lpthread','-lm','-ldl']
 link=[str(mpi/'mpiicc'),*common,*[str(out/'objects'/u['object'])for u in plan['units']],*libs,'-o',str(out/'openmx-cpu')]
 return commands,link

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--selection',type=Path)
 p.add_argument('--profile',choices=['stock','soc','unified'],default='unified');p.add_argument('--jobs',type=int,default=8)
 p.add_argument('--oneapi-prefix',type=Path,default=Path(DEFAULT_ONEAPI));p.add_argument('--cpu-target',default='CORE-AVX2',help='ICC target, default portable Intel AVX2; HOST is native-only')
 p.add_argument('--execute',action='store_true');a=p.parse_args()
 if not 1<=a.jobs<=16:p.error('--jobs must be 1..16')
 selection=a.selection
 if a.profile=='unified'and selection is None:selection=REPO/'cpu/overrides/nsoc/UNIFIED_SELECTION.json'
 try:sel=json.loads(selection.read_text())if selection else None;plan=source_plan(sel, a.profile)
 except(ValueError,OSError)as e:p.error(str(e))
 out=a.output.resolve()
 if out.exists():p.error('output must not exist; every build compiles into an empty directory')
 commands,link=command_plan(plan,out,a.oneapi_prefix.resolve(),a.cpu_target)
 inputs={u['source']:sha(local(u['source']))for u in plan['units']}
 # Include all local headers/templates in declared source closure.
 for d in plan['include_dirs']:
  for f in local(d).rglob('*'):
   if f.is_file()and f.suffix.lower()in ['.c','.h','.f','.inc','.f90']:inputs[str(f.relative_to(REPO))]=sha(f)
 runtime_sources={str(f.relative_to(REPO)):sha(f)for f in runtime_resource_files()}if plan['capabilities'].get('input_symmetry')else{}
 manifest={'schema':1,'state':'planned','hostname':socket.gethostname(),'cwd':str(REPO),'environment_prefix':str(a.oneapi_prefix.resolve()),'python_prefix':os.sys.prefix,'builder_sha256':sha(__file__),'cpu_target':a.cpu_target,'profile':a.profile,'selection':str(selection)if selection else None,'source_manifest':str(BASE.relative_to(REPO)),'sources_sha256':inputs,'capabilities':plan['capabilities'],'external_openmx_objects':[],'external_elpa_objects':[],'build_environment':{'I_MPI_ROOT':str(a.oneapi_prefix.resolve()/'mpi/2021.7.1')},'jobs':a.jobs,'compile_commands':commands,'link_command':link,'binary':{'path':str(out/'openmx-cpu'),'sha256':None},'numerics_validated':False}
 manifest['runtime_resources_sha256']=runtime_sources
 manifest['runtime_python']=sys.executable
 if not a.execute:print(json.dumps(manifest,indent=2));return
 if socket.gethostname().split('.')[0].startswith('login'):p.error('compile on an authorized compute node')
 (out/'objects').mkdir(parents=True);(out/'logs').mkdir()
 if runtime_sources:copy_runtime_resources(out)
 shell='source '+shlex.quote(str(a.oneapi_prefix.resolve()/'setvars.sh'))+' >/dev/null 2>&1; export I_MPI_ROOT='+shlex.quote(str(a.oneapi_prefix.resolve()/'mpi/2021.7.1'))+'; exec "$@"'
 def run(cmd,label):
  start=time.monotonic();q=subprocess.run(['/bin/bash','-c',shell,'cpu-build',*cmd],cwd=out/'objects',stdin=subprocess.DEVNULL,capture_output=True,text=True)
  (out/'logs'/(label+'.stdout.log')).write_text(q.stdout);(out/'logs'/(label+'.stderr.log')).write_text(q.stderr)
  rec={'label':label,'command':cmd,'returncode':q.returncode,'seconds':time.monotonic()-start}
  if q.returncode:raise RuntimeError(label+' failed; '+str(out/'logs'/(label+'.stderr.log')))
  return rec
 start=time.monotonic();manifest['state']='building';(out/'build-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');records=[]
 try:
  # F90 modules are built in stock makefile dependency order; C translation units are independent.
  with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs)as pool:
   futures=[pool.submit(run,c['command'],c['label'])for c in commands if c['language']=='c']
   for c in commands:
    if c['language']=='fortran':records.append(run(c['command'],c['label']))
   for future in concurrent.futures.as_completed(futures):records.append(future.result())
  records.append(run(link,'link'))
  tools=[str(a.oneapi_prefix/'mpi/2021.7.1/bin/mpiicc'),str(a.oneapi_prefix/'mpi/2021.7.1/bin/mpiifort'),str(a.oneapi_prefix/'mpi/2021.7.1/bin/mpiexec')]
  manifest['compiler_versions']={t:run([t,'--version'],'version-'+Path(t).name)for t in tools}
  q=subprocess.run(['/bin/bash','-c',shell,'cpu-build','ldd',str(out/'openmx-cpu')],stdin=subprocess.DEVNULL,capture_output=True,text=True);(out/'logs/ldd.txt').write_text(q.stdout+q.stderr)
  if q.returncode or 'not found'in q.stdout:raise RuntimeError('linked runtime unresolved; see logs/ldd.txt')
  manifest['loaded_runtime_abi']=q.stdout.splitlines()
  manifest['compiler_version_text']={Path(t).name:(out/'logs'/('version-'+Path(t).name+'.stdout.log')).read_text().strip()for t in tools}
  manifest['binary']['sha256']=sha(out/'openmx-cpu');manifest['state']='built_pending_native_validation'
 except Exception as e:manifest.update(state='failed',error=str(e));raise
 finally:
  manifest.update(records=records,seconds=time.monotonic()-start)
  (out/'build-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({'state':manifest['state'],'binary':manifest['binary'],'seconds':manifest['seconds']},indent=2))
if __name__=='__main__':main()
