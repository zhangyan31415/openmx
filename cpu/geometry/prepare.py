#!/usr/bin/env python3
"""Prepare one explicit closed subgroup for both native NSOC and one-step SOC.

No external tapw/OpenMX source tree, geometry modification, or solver launch.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
TOOLS=HERE/'metadata'
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
 if soc_route not in ('generic','legacy','auto'):raise ValueError('SOC route must be generic, legacy or auto')
 if export_policy not in ('raw','symmetrize'):raise ValueError('SOC export policy must be raw or symmetrize')
 if soc_route=='legacy'and group not in ('auto','D3','C3'):
  raise ValueError('legacy SOC route supports only C3 or D3; use generic for '+str(group))

def automatic_soc_route(group,structure,tolerance,operations=None):
 """Use the verified slab solver only in its supported angular convention."""
 if group not in ('C3','D3'):return 'generic'
 q=np.array([[-.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,-.5,0],[0,0,1.]])
 cell=np.asarray(structure['cell']);bound=float(tolerance)*1.25
 flat=abs(cell[0,2])+abs(cell[1,2])+abs(cell[2,0])+abs(cell[2,1])
 if flat>4*bound or np.linalg.norm(q@cell[0]-cell[1])>bound:return 'generic'
 if operations is not None:
  rotations=[np.asarray(g['R_cart_orthogonal'])for g in operations]
  if not any(np.max(abs(r-q))<1e-3 for r in rotations):return 'generic'
  if group=='D3'and not any(np.linalg.det(r)>.99 and r[2,2]<-.99 for r in rotations):return 'generic'
 return 'legacy'

def prepare_soc_artifacts(soc_input,selected_path,out,kind,spinor,neighbor_settings,*,soc_route='generic',export_policy='raw',realspace_map=None):
 """Prepare exactly one SOC representation and declare native export semantics.

 This returns environment additions, not a complete runtime profile. The
 caller must clear inherited OpenMX selectors before applying that profile.
 """
 _soc_options(kind,soc_route,export_policy)
 out=Path(out)
 if soc_route=='generic':
  sbasis=out/'spinor.gsoc';smeta=out/'spinor.gsoc.json'
  spinor_info=spinor.build(soc_input,selected_path,sbasis,smeta,kind)
  files={'soc_basis':sbasis,'soc_metadata':smeta}
  env={'OPENMX_CPU_GSOC':'1','OPENMX_CPU_GSOC_BASIS':str(sbasis),**neighbor_settings}
  export_HS='raw'
  if export_policy=='symmetrize':
   helper=module('unified_generic_soc_export',HERE/'generic_soc_export.py')
   rmap=out/'operations.socrs';rmeta=out/'operations.socrs.json'
   helper.build(soc_input,spinor_info,rmap,rmeta)
   files.update(soc_native_map=rmap,soc_native_metadata=rmeta)
   env.update(OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE='1',OPENMX_CPU_SOC_NATIVE_MAP=str(rmap),OPENMX_CPU_GSOC_PROJECT='1')
   export_HS='symmetrized'
 else:
  if kind not in ('D3','C3'):raise ValueError('selected subgroup is unsupported by legacy SOC')
  legacy=module('unified_legacy_soc_prepare',HERE/'legacy_soc.py')
  info=legacy.prepare(soc_input,selected_path,out/'legacy_soc',kind,export_policy=export_policy)
  files={'soc_basis':Path(info['basis']),'soc_metadata':out/'legacy_soc/LEGACY_SOC.json',
         'soc_native_map':Path(info['map'])}
  env={**info['env_additions'],**neighbor_settings};export_HS=info['export_HS']
 return {'route':soc_route,'export_policy':export_policy,'export_HS':export_HS,
         'files':files,'env_additions':env}

def tolerance_ladder(limit, supplied=None):
 """Include the user's exact upper bound, without expanding the solve scope."""
 import math
 limit=float(limit)
 if not math.isfinite(limit) or not 0<limit<=.02:
  raise ValueError('geometry tolerance must be finite in (0,0.02] Angstrom')
 if supplied is None:
  values=[1e-6,3e-6,1e-5,2e-5,3e-5,5e-5,1e-4,
          2e-4,3e-4,5e-4,1e-3,2e-3,3e-3,5e-3,1e-2,2e-2]
 else:
  values=[float(v)for v in supplied]
  if not values or any(not math.isfinite(v)or v<=0 for v in values)or any(a>=b for a,b in zip(values,values[1:])):
   raise ValueError('custom tolerance ladder must be finite, positive and strictly increasing')
 values=[v for v in values if v<limit]
 # A second upper-range probe allows a new symmetry plateau to be observed.
 if not values or values[-1]<limit/2 or limit>1e-4:values.append(limit*.75)
 values.append(limit)
 return sorted(set(values))


def failure_diagnostic(input_path, full, ladder, limit, detector, spinor, group):
 """Inspect looser candidates, but never authorize them for a native solve."""
 import math
 def brief(item):
  pg=item.get('point_group',{})
  return {'symprec_Ang':item.get('symprec_Ang'), 'valid':bool(item.get('valid')),
          'point_group':pg.get('symbol')if isinstance(pg,dict)else pg,
          'operation_count':item.get('operation_count',0),
          'max_geometry_residual_Ang':(float(item['max_geometry_residual_Ang'])if item.get('max_geometry_residual_Ang')is not None and math.isfinite(item['max_geometry_residual_Ang'])else None),
          'error':item.get('error',item.get('reason'))}
 scan=list(full.get('scan')or[full.get('selected',{})])
 report={'schema':1,'requested_group':group,'requested_max_Ang':limit,
         'solve_scan_Ang':list(ladder),'solve_scan':[brief(r)for r in scan],
         'diagnostic_scan':[],'candidate':None,'absence_proven':False,
         'runtime_launch_allowed':False,'geometry_modified':False,
         'tolerance_automatically_increased':False}
 nontrivial=[r for r in scan if r.get('valid')and r.get('operation_count',0)>1]
 if nontrivial:
  report['classification']='unsupported_subgroup'
  for record in nontrivial:
   try:kind,ops=select_subgroup(record,group,spinor)
   except ValueError:continue
   report['classification']='unstable_or_unusable_selection'
   break
  return report
 report['classification']='operation_validation_failed'if any(not r.get('valid')for r in scan)else'identity_only_in_scanned_range'
 probes=sorted({min(.02,limit*factor)for factor in (3,10,30)if min(.02,limit*factor)>limit})
 structure=detector.parse_openmx(input_path)
 candidates=[]
 for tol in probes:
  try:
   record=detector.detect_at_tolerance(structure,tol,1.25)
   item=brief(record)
   residual=record.get('max_geometry_residual_Ang')
   if record.get('valid')and record.get('operation_count',0)>1 and residual is not None and math.isfinite(residual)and residual<=tol:
    try:kind,ops=select_subgroup(record,group,spinor)
    except ValueError:pass
    else:
     item['supported_subgroup']=kind
     candidates.append({'symprec_Ang':tol,'group':kind,'operation_count':len(ops),
                        'max_geometry_residual_Ang':residual})
  except (ValueError,RuntimeError)as exc:
   item={'symprec_Ang':tol,'valid':False,'error':str(exc)}
  report['diagnostic_scan'].append(item)
 if candidates:
  report['classification']='larger_tolerance_candidate'
  report['candidate']=min(candidates,key=lambda r:(r['symprec_Ang'],-r['operation_count']))
 return report


def diagnostic_message(report):
 limit=report['requested_max_Ang'];kind=report['classification']
 lines=[f'No usable symmetry subgroup selected within scf.Symmetry.Tolerance={limit:.6e} Ang; diagnostic={kind}.']
 if kind=='larger_tolerance_candidate':
  c=report['candidate']
  lines.append(f'Diagnostic-only scan found a validated {c["group"]} candidate at {c["symprec_Ang"]:.6e} Ang (maximum position residual {c["max_geometry_residual_Ang"]:.6e} Ang).')
  lines.append(f'If that geometric approximation is appropriate, set scf.Symmetry.Tolerance {c["symprec_Ang"]:.6e} and rerun. A candidate alone does not establish physical accuracy or a stable tolerance platform.')
 elif kind=='identity_only_in_scanned_range':
  lines.append('Only the identity was validated in the solve scan; this does not prove that the structure has no symmetry. Inspect the geometry, or use scf.Symmetry off for an unrestricted calculation.')
 elif kind=='operation_validation_failed':
  lines.append('An operation/mapping validation failed; this is not a conclusion that the structure has no symmetry. Inspect the diagnostic mapping errors and geometry before increasing the tolerance.')
 elif kind=='unsupported_subgroup':
  lines.append('Nontrivial operations were detected, but no supported subgroup matches the request. Increasing the tolerance is not necessarily a remedy; inspect the group or use scf.Symmetry off.')
 else:
  lines.append('Nontrivial candidates occur in the scan but the selected operation set is unstable or unusable. Inspect the tolerance platform and mapping validation.')
 if report['diagnostic_scan']:
  tried=', '.join(f'{r["symprec_Ang"]:.6e}'for r in report['diagnostic_scan'])
  lines.append('Additional diagnostic-only tolerances (Ang): '+tried+'.')
 lines.append('Native SCF has not started. The geometry and requested solve tolerance were not changed automatically.')
 return '\n'.join(lines)


def stop_with_diagnostic(input_path,out,full,ladder,limit,detector,spinor,group,cause):
 report=failure_diagnostic(input_path,full,ladder,limit,detector,spinor,group)
 report['original_error']=str(cause)
 path=out/'SYMMETRY_DIAGNOSTIC.json';dump(path,report)
 message=diagnostic_message(report)+'\nOriginal validation error: '+str(cause)+'\nDiagnostic details: '+str(path)
 (out/'symmetry-diagnostic.txt').write_text(message+'\n')
 raise ValueError(message)from cause


def prepare(nsoc_input,soc_input,output,group='auto',tolerances=None,atom_tol=5e-5,project_nsoc=False,*,soc_route='generic',export_policy='raw',detection_mode='auto',operations_file=None):
 _soc_options(group,soc_route,export_policy)
 nsoc_input=Path(nsoc_input).resolve();soc_input=Path(soc_input).resolve();out=Path(output).resolve()
 workflow=module('prepare_workflow_contract',REPO/'cpu/tools/workflow.py')
 workflow.validate_inputs(nsoc_input.read_text(),soc_input.read_text())
 if out.exists():raise ValueError('output must be a new metadata directory')
 out.mkdir(parents=True)
 detector=module('unified_detector',TOOLS/'symmetry_detect.py')
 spinor=module('unified_spinor_prepare',HERE/'prepare_soc_basis.py')
 realspace=module('unified_realspace_prepare',TOOLS/'prepare_realspace_d3.py')
 try:
  if detection_mode not in ('auto','user'):raise ValueError('detection mode must be auto or user')
  ladder=tolerance_ladder(atom_tol,tolerances)
  full=detector.run(nsoc_input,out/'detected.json',ladder,5e-6,1.25,
    mode=detection_mode,requested_symprec=atom_tol if detection_mode=='user'else None,
    operations_file=Path(operations_file)if operations_file else None)
  checked=full.get('selected',{})
  try:
   if not checked.get('valid'):raise ValueError('no valid detected operation set')
   kind,operations=select_subgroup(checked,group,spinor)
  except ValueError as exc:
   stop_with_diagnostic(nsoc_input,out,full,ladder,atom_tol,detector,spinor,group,exc)
  if operations_file and len(operations)!=checked['operation_count']:
   raise ValueError('OperationsFile must be exactly the requested supported subgroup; refusing silent subgroup reduction')
  structure=detector.parse_openmx(nsoc_input)
  symprec=float(checked['symprec_Ang'])
  selected=detector._build_operation_set(structure,[g['W']for g in operations],[g['t']for g in operations],
    symprec,1.25,source='unified_explicit_'+kind)
  if not selected.get('valid')or selected['max_geometry_residual_Ang']>atom_tol:
   stop_with_diagnostic(nsoc_input,out,full,ladder,atom_tol,detector,spinor,group,ValueError('selected subgroup fails fixed-geometry validation'))
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
  chosen_route=automatic_soc_route(kind,structure,symprec,selected['operations'])if soc_route=='auto'else soc_route
  soc=prepare_soc_artifacts(soc_input,selected_path,out,kind,spinor,neighbor['settings'],
                            soc_route=chosen_route,export_policy=export_policy)
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
    'geometry_modified':False,'solve_tolerance_ladder_Ang':ladder,'requested_geometry_tolerance_Ang':atom_tol,'files':{k:str(v)for k,v in files.items()},'sha256':{k:sha(v)for k,v in files.items()},
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
