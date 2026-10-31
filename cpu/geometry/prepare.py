#!/usr/bin/env python3
"""Prepare one explicit closed subgroup for both native NSOC and one-step SOC.

No external tapw/OpenMX source tree, geometry modification, or solver launch.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
TOOLS=REPO/'local-records/cpu_d3/historical/geometry_tools'
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):
 def convert(x):
  if isinstance(x,np.ndarray):return x.tolist()
  if isinstance(x,np.generic):return x.item()
  if isinstance(x,Path):return str(x)
  raise TypeError(type(x).__name__)
 Path(path).write_text(json.dumps(value,indent=2,default=convert)+'\n')

def select_subgroup(checked,requested,spinor):
 if requested in ('D2','C2'):
  kind,ops,_,_=spinor._subgroup(checked,requested)
 elif requested in ('Ci','Cs'):
  kind,ops,_,_=spinor._improper_subgroup(checked,requested)
 else:
  try:kind,ops,_,_=spinor._cyclic_subgroup(checked,requested)
  except ValueError:
   if requested!='auto':raise
   kind,ops,_,_=spinor._improper_subgroup(checked,'auto')
 return kind,ops

def _soc_options(group,soc_route,export_policy):
 if soc_route not in ('generic','legacy'):raise ValueError('SOC route must be generic or legacy')
 if export_policy not in ('raw','symmetrize'):raise ValueError('SOC export policy must be raw or symmetrize')
 if soc_route=='generic'and export_policy!='raw':
  raise ValueError('generic SOC currently exports raw H/S; symmetrized export requires explicit legacy C3/D3')
 if soc_route=='legacy'and group not in ('auto','D3','C3'):
  raise ValueError('legacy SOC route supports only C3 or D3; use generic for '+str(group))

def prepare_soc_artifacts(soc_input,selected_path,out,kind,spinor,neighbor_settings,*,soc_route='generic',export_policy='raw'):
 """Prepare exactly one SOC representation and declare native export semantics.

 This returns environment additions, not a complete runtime profile. The
 caller must clear inherited OpenMX selectors before applying that profile.
 """
 _soc_options(kind,soc_route,export_policy)
 out=Path(out)
 if soc_route=='generic':
  sbasis=out/'spinor.gsoc';smeta=out/'spinor.gsoc.json'
  spinor.build(soc_input,selected_path,sbasis,smeta,kind)
  files={'soc_basis':sbasis,'soc_metadata':smeta}
  env={'OPENMX_CPU_GSOC':'1','OPENMX_CPU_GSOC_BASIS':str(sbasis),**neighbor_settings}
  export_HS='raw'
 else:
  if kind not in ('D3','C3'):raise ValueError('selected subgroup is unsupported by legacy SOC')
  legacy=module('unified_legacy_soc_prepare',HERE/'legacy_soc.py')
  info=legacy.prepare(soc_input,selected_path,out/'legacy_soc',kind,export_policy=export_policy)
  files={'soc_basis':Path(info['basis']),'soc_metadata':out/'legacy_soc/LEGACY_SOC.json',
         'soc_native_map':Path(info['map'])}
  env={**info['env_additions'],**neighbor_settings};export_HS=info['export_HS']
 return {'route':soc_route,'export_policy':export_policy,'export_HS':export_HS,
         'files':files,'env_additions':env}

def prepare(nsoc_input,soc_input,output,group='auto',tolerances=None,atom_tol=5e-5,project_nsoc=False,*,soc_route='generic',export_policy='raw'):
 _soc_options(group,soc_route,export_policy)
 nsoc_input=Path(nsoc_input).resolve();soc_input=Path(soc_input).resolve();out=Path(output).resolve()
 workflow=module('prepare_workflow_contract',REPO/'cpu/tools/workflow.py')
 workflow.validate_inputs(nsoc_input.read_text(),soc_input.read_text())
 if out.exists():raise ValueError('output must be a new metadata directory')
 out.mkdir(parents=True)
 detector=module('unified_detector',TOOLS/'symmetry_detect.py')
 spinor=module('unified_spinor_prepare',REPO/'local-records/cpu_d3/generic_soc/tools/prepare_soc_basis.py')
 realspace=module('unified_realspace_prepare',TOOLS/'prepare_realspace_d3.py')
 ladder=tolerances or [1e-6,3e-6,1e-5,2e-5,3e-5,5e-5,1e-4]
 try:
  full=detector.run(nsoc_input,out/'detected.json',ladder,5e-6,1.25)
  checked=full.get('selected',{})
  if not checked.get('valid'):raise ValueError('no valid detected operation set')
  kind,operations=select_subgroup(checked,group,spinor)
  structure=detector.parse_openmx(nsoc_input)
  symprec=float(checked['symprec_Ang'])
  selected=detector._build_operation_set(structure,[g['W']for g in operations],[g['t']for g in operations],
    symprec,1.25,source='unified_explicit_'+kind)
  if not selected.get('valid')or selected['max_geometry_residual_Ang']>atom_tol:
   raise ValueError('selected subgroup fails fixed-geometry validation')
  metadata={'input':str(nsoc_input),'input_sha256':sha(nsoc_input),'selected':selected,'selected_group':kind,
            'source_metadata':str(out/'detected.json'),'source_metadata_sha256':sha(out/'detected.json'),
            'geometry_modified':False,'natoms':structure['natoms'],'species':structure['species'],'species_basis':structure['species_basis']}
  selected_path=out/'selected.symmetry.json';dump(selected_path,metadata)
  nsbasis=out/'scalar.nsym';nsmeta=out/'scalar.nsym.json'
  if kind=='D3':
   scalar=module('unified_paired_d3',REPO/'cpu/overrides/nsoc/metadata/prepare_paired_d3.py')
   scalar.build(nsoc_input,selected_path,nsbasis,nsmeta);version=3
  elif kind=='C3':
   scalar=module('unified_paired_c3',REPO/'cpu/overrides/nsoc/metadata/prepare_paired_c3.py')
   scalar.build(nsoc_input,selected_path,nsbasis,nsmeta);version=4
  else:
   scalar=module('unified_generic_scalar',TOOLS/'prepare_nsoc_auto.py')
   scalar.build(nsoc_input,selected_path,nsbasis,nsmeta);version=2
  rmap=out/'operations.rsd3';rmeta=out/'operations.rsd3.json'
  realspace.generate(nsoc_input,rmap,provenance_path=rmeta,symprec=symprec,atom_tol=atom_tol,
                     include_improper=True,operations_metadata=metadata)
  support=module('unified_neighbor_support',TOOLS/'neighbor_support.py')
  neighbor_path=out/'neighbor_support.json'
  neighbor=support.prepare(nsoc_input,rmeta,neighbor_path,max_padding=.01,runtime_cwd=nsoc_input.parent)
  soc=prepare_soc_artifacts(soc_input,selected_path,out,kind,spinor,neighbor['settings'],
                            soc_route=soc_route,export_policy=export_policy)
  if project_nsoc and kind not in ('D3','C3'):
   raise ValueError('generic real-space H/S projection is not production-validated; explicit D3/C3 only')
  adds={'OPENMX_CPU_NSOC_SYM':'1','OPENMX_CPU_NSOC_FULL_PROJECT':'1','OPENMX_CPU_NSOC_BASIS':str(nsbasis),
   'OPENMX_CPU_NSOC_GENERIC_COV_STRICT':'1','OPENMX_CPU_REALSPACE_D3':'scf'if project_nsoc else'off',
   'OPENMX_CPU_REALSPACE_MAP':str(rmap),'OPENMX_CPU_REALSPACE_CHECK':'0'}
  adds.update(neighbor['settings'])
  files={'selected_symmetry':selected_path,'nsoc_basis':nsbasis,'nsoc_metadata':nsmeta,
         **soc['files'],'realspace_map':rmap,'realspace_metadata':rmeta,'neighbor_support':neighbor_path}
  result={'schema':1,'status':'prepared','requested_group':group,'selected_group':kind,
    'detected_point_group':checked['point_group'],'NSOC_basis_version':version,'NSOC_project_explicit':project_nsoc,
    'NSOC_realspace_HS_projection_explicit':project_nsoc,
    'NSOC_export_HS':'PROJECTED_REALSPACE_H_S'if project_nsoc else'RAW_REALSPACE_H_S',
    'NSOC_density_export_projection':False,
    'SOC_route':soc['route'],'SOC_export_policy':soc['export_policy'],'SOC_export_HS':soc['export_HS'],
    'scope':'fixed nonmagnetic Gamma NSOC -> one SOC; metadata only, no runtime covariance or physical accuracy claim',
    'geometry_modified':False,'files':{k:str(v)for k,v in files.items()},'sha256':{k:sha(v)for k,v in files.items()},
    'NSOC_environment_additions':adds,'SOC_environment_additions':soc['env_additions'],
    'full_runtime_profile_required':True,'runtime_HS_covariance_checked':False}
  dump(out/'PREPARED.json',result);return result
 except Exception as exc:
  dump(out/'PREPARED.json',{'schema':1,'status':'failed','error':str(exc),'runtime_launch_allowed':False})
  raise

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--nsoc-input',type=Path,required=True);p.add_argument('--soc-input',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--group',choices=['auto','D3','C3','D2','C2','D4','C4','D6','C6','Ci','Cs'],default='auto')
 p.add_argument('--symprec-ladder',type=float,nargs='+');p.add_argument('--geometry-tol',type=float,default=5e-5)
 p.add_argument('--project-nsoc',action='store_true',help='project real-space H/S during NSOC and export projected H/S (C3/D3 only); DM/EDM hooks are read-only checks')
 p.add_argument('--soc-route',choices=['generic','legacy'],default='generic')
 p.add_argument('--soc-export-policy','--export-policy',dest='soc_export_policy',choices=['raw','symmetrize'],default='raw')
 a=p.parse_args()
 try:r=prepare(a.nsoc_input,a.soc_input,a.output,a.group,a.symprec_ladder,a.geometry_tol,a.project_nsoc,
               soc_route=a.soc_route,export_policy=a.soc_export_policy)
 except (ValueError,OSError,RuntimeError)as e:p.error(str(e))
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
