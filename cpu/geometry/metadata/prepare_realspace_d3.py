#!/usr/bin/env python3
"""Prepare a closed scalar point-group operation set for NSOC projection.

Binary v1 is little endian, with no padding:
 int32[8]: magic=0x52534433, version=1, natoms,N,operation_count,nspecies,0,0;
 float64 cell_Ang[3,3], xyz_Ang[natoms,3];
 int32 atom_species[natoms], species_norb[nspecies]; then each operation:
 int32 W[3,3], perm[natoms], image_shift[natoms,3];
 float64 D_species[norb,norb] for each species (all arrays row major).

Geometry is parsed without alteration by the pinned release package. The
six-operation D3 case keeps the established exact Q3/Q2 angular convention;
other scalar groups use the package orbital blocks after orthogonality and
composition checks. Integer atom/image mappings remain package values. There
are no spin factors. Neighbor keys use ((Rx,Ry,Rz), source_atom, neighbor_atom),
all zero based. This utility's small Python projector is a mathematical
reference, not the native execution.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import socket
import sys
import time

_MODULE_IMPORT_STARTED=time.perf_counter()
import numpy as np

from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
PACKAGE_ROOT=HERE.parent
NOMINAL_HELPER=HERE/'prepare_nsoc_d3.py'
_spec=importlib.util.spec_from_file_location('native_geometry_detector',HERE/'symmetry_detect.py')
_detector=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_detector)
ANGULAR_MODULE=PACKAGE_ROOT/'openmx_orbitals.py'
MODULE_IMPORT_SECONDS=time.perf_counter()-_MODULE_IMPORT_STARTED


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def transform_key(key, operation):
    """Active periodic action; preserve the source-to-target image convention."""
    r,i,j=key
    p=operation['perm']; image=operation['image_shift']
    new_r=operation['W']@np.asarray(r,dtype=np.int64)+image[j]-image[i]
    return (tuple(int(v) for v in new_r),int(p[i]),int(p[j]))


def validate_closure(blocks, operations, *, require_reverse=True):
    """Reject missing symmetry or Hermitian reverse support; return key count.

    Values are deliberately ignored: a present all-zero block remains part of
    finite-range support. Missing blocks are never silently supplied as zeros.
    """
    keys=set(blocks)
    for key in sorted(keys):
        if require_reverse:
            r,i,j=key
            reverse=(tuple(-int(v) for v in r),j,i)
            if reverse not in keys:
                raise ValueError(f'missing Hermitian reverse closure: {key} -> {reverse}')
        for g,operation in enumerate(operations):
            target=transform_key(key,operation)
            if target not in keys:
                raise ValueError(f'missing D3 neighbor closure: g={g}, {key} -> {target}')
    return len(keys)


def transform_blocks(blocks, operation, atom_species):
    """Independent small-array reference for real NSOC H and S blocks."""
    out={}
    for key,block in blocks.items():
        _,i,j=key
        left=operation['rotations'][int(atom_species[i])]
        right=operation['rotations'][int(atom_species[j])]
        target=transform_key(key,operation)
        if target in out:
            raise ValueError(f'non-bijective block action at {target}')
        out[target]=left@block@right.T
    return out


def project_blocks(blocks, operations, atom_species):
    """Average a closed scalar operation set on complete finite-range support."""
    if not operations or len(operations)>48:
        raise ValueError('the operation set must contain 1 through 48 operations')
    validate_closure(blocks,operations)
    out={key:np.zeros_like(value,dtype=np.result_type(value,np.float64)) for key,value in blocks.items()}
    for operation in operations:
        for key,block in transform_blocks(blocks,operation,atom_species).items():
            out[key]+=block/len(operations)
    return out


def group_composition_error(operations, *, orbital_tol=1e-12):
    """Verify integer atom/image action and noncommuting orbital representation."""
    if not operations or len(operations)>48:
        raise ValueError('the operation set must contain 1 through 48 operations')
    n=len(operations[0]['perm']); identity=np.arange(n)
    if len({tuple(g['W'].ravel()) for g in operations})!=len(operations):
        raise ValueError('duplicate fractional rotations')
    largest=0.
    for g in operations:
        if not np.array_equal(np.sort(g['perm']),identity):
            raise ValueError('atom mapping is not a permutation')
        for h in operations:
            w=g['W']@h['W']
            match=[k for k in operations if np.array_equal(k['W'],w)]
            if len(match)!=1:
                raise ValueError('point-group rotations do not close')
            k=match[0]
            if not np.array_equal(g['perm'][h['perm']],k['perm']):
                raise ValueError('point-group atom permutations do not compose')
            shifts=h['image_shift']@g['W'].T+g['image_shift'][h['perm']]
            # Spglib can choose group representatives differing by a common
            # lattice translation. It cancels exactly from image_j-image_i.
            difference=shifts-k['image_shift']
            if not np.all(difference==difference[0]):
                raise ValueError('point-group periodic relative image shifts do not compose')
            for dg,dh,dk in zip(g['rotations'],h['rotations'],k['rotations']):
                largest=max(largest,float(np.max(abs(dg@dh-dk))))
    if largest>orbital_tol:
        raise ValueError(f'orbital point-group composition error {largest}')
    return largest


def nominal_rotations(input_path):
    """Use the established exact Q3/Q2 without changing source geometry."""
    spec=importlib.util.spec_from_file_location('established_nsoc_d3_geometry',NOMINAL_HELPER)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    _,_,_,q3,_,_,q2,_,_=helper.geometry(Path(input_path).read_text())
    return [np.eye(3),q3,q3@q3,q2,q3@q2,q3@q3@q2]


def prepare(input_path, *, symprec=1e-4, atom_tol=5e-6, include_improper=False, operations_metadata=None):
    """Use exactly the supplied validated operation set; never redetect a larger group."""
    raw=_detector.parse_openmx(Path(input_path))
    if raw['spin_mode']!='none':
        raise ValueError('real-space metadata requires nonmagnetic scalar NSOC')
    explicit=operations_metadata is not None
    if isinstance(operations_metadata,(str,Path)):
        operations_metadata=json.loads(Path(operations_metadata).read_text())
    if explicit:
        selected=operations_metadata.get('selected',operations_metadata)
        entries=selected.get('operations',[])
        if not selected.get('valid') or not entries:
            raise ValueError('a validated nonempty operation set is required')
    else:
        selected=_detector.detect_at_tolerance(raw,symprec,max(1.25,atom_tol/symprec))
        entries=selected.get('operations',[])
        if not include_improper:
            entries=[op for op in entries if np.linalg.det(op['W'])>.5]
    checked=_detector._build_operation_set(raw,[op['W']for op in entries],
        [op.get('t',op.get('translation_frac',[0,0,0]))for op in entries],
        symprec,max(1.25,atom_tol/symprec),source='explicit_metadata'if explicit else'repository_detector')
    if not checked.get('valid'):
        raise ValueError('operation set failed geometry/closure/AO validation: '+str(checked.get('error',checked.get('closure'))))
    entries=checked['operations']
    if not entries or len(entries)>48:
        raise ValueError('point-group operation count must be1..48')
    entries.sort(key=lambda op:0 if np.array_equal(op['W'],np.eye(3,dtype=int)) and np.max(abs(np.asarray(op['t'])-np.rint(op['t'])))<1e-10 else 1)
    species=list(raw['species_basis'])
    counts={name:sum((2*'spdf'.index(shell)+1)*num for shell,num in _detector._basis_shells(basis))
            for name,basis in raw['species_basis'].items()}
    atom_species=np.array([species.index(name)for name in raw['species']],dtype=np.int32)
    species_norb=np.array([counts[name]for name in species],dtype=np.int32)
    structure=SimpleNamespace(lattice=raw['cell'],species_basis=raw['species_basis'],spin_mode='none',spinful=False,
        nwann=sum(counts[name]for name in raw['species']),
        atoms=[SimpleNamespace(species=name,frac=raw['frac'][i],cart=raw['cart'][i])for i,name in enumerate(raw['species'])])
    is_d3=len(entries)==6 and sorted(round(float(np.trace(op['R_cart_orthogonal'])))for op in entries)==[-1,-1,-1,0,0,3]
    nominal=None
    if is_d3:
        try:nominal=nominal_rotations(input_path)
        except ValueError:pass # Arbitrary orientation uses the validated common AO representation.
    if nominal is not None:
        # A hexagonal cell admits two distinct D3 C2-axis families. A nominal
        # gauge is valid only for the actually selected operation set, not
        # merely because the cell admits some conventional D3 orientation.
        matches=[int(np.argmin([np.max(abs(v-np.asarray(op['R_cart_orthogonal'])))for v in nominal]))for op in entries]
        if len(set(matches))!=len(entries)or any(np.max(abs(nominal[k]-np.asarray(op['R_cart_orthogonal'])))>1e-6 for k,op in zip(matches,entries)):
            nominal=None
    operations=[];transports=[];used=set()
    for index,entry in enumerate(entries):
        w=np.asarray(entry['W'],dtype=int);t=np.asarray(entry['t'],dtype=float)
        q=np.asarray(entry['R_cart_orthogonal'],dtype=float);match=-1;qnom=q
        if nominal is not None:
            distances=[float(np.max(abs(v-q)))for v in nominal];match=int(np.argmin(distances))
            if distances[match]>1e-6 or match in used:raise ValueError('nominal D3 operation mismatch')
            used.add(match);qnom=nominal[match]
        rotations=[];raw_rotations=[]
        for name in species:
            basis=raw['species_basis'][name]
            a=_detector._orbital_rotation_block(basis,q)
            d=_detector._orbital_rotation_block(basis,qnom)
            if max(np.max(abs(a.imag),initial=0),np.max(abs(d.imag),initial=0))>1e-12:
                raise ValueError('scalar AO rotation must be real')
            if np.max(abs(d.real.T@d.real-np.eye(len(d))),initial=0)>1e-12:
                raise ValueError('AO rotation is not orthogonal')
            rotations.append(d.real);raw_rotations.append(a.real)
        perm=np.asarray(entry['atom_permutation'],dtype=np.int32)
        shifts=np.asarray(entry['image_shift'],dtype=np.int32)
        operation=SimpleNamespace(index=index,rotation_frac=w,translation_frac=t,rotation_cart=q,
                                  translation_cart=raw['cell'].T@t)
        transports.append(SimpleNamespace(operation=operation,atom_mappings=perm,image_shifts=shifts))
        operations.append(dict(W=w,perm=perm,image_shift=shifts,rotations=rotations,
                               rotation_cart_nominal=qnom,rotation_cart_cell=q,
                               cell_derived_rotations=raw_rotations,nominal_index=match))
    err=group_composition_error(operations)
    for op in operations:
        if not np.array_equal(atom_species[op['perm']],atom_species):raise ValueError('operation changes species')
    return structure,species,atom_species,species_norb,transports,operations,err


def generate(input_path, output_path, *, provenance_path=None, symprec=1e-4, atom_tol=5e-6, include_improper=False, operations_metadata=None):
    started=time.perf_counter()
    input_path=Path(input_path).resolve(); output_path=Path(output_path).resolve()
    structure,species,atom_species,species_norb,transports,operations,err=prepare(input_path,symprec=symprec,atom_tol=atom_tol,include_improper=include_improper,operations_metadata=operations_metadata)
    cell=np.asarray(structure.lattice); xyz=np.array([a.cart for a in structure.atoms])
    chunks=[np.array([0x52534433,1,len(xyz),structure.nwann,len(operations),len(species),0,0],dtype='<i4').tobytes(),
            cell.astype('<f8').tobytes(),xyz.astype('<f8').tobytes(),atom_species.astype('<i4').tobytes(),species_norb.astype('<i4').tobytes()]
    operation_payload=[]
    for op,t in zip(operations,transports):
        for key in ('W','perm','image_shift'):chunks.append(np.asarray(op[key],dtype='<i4').tobytes())
        chunks.extend(d.astype('<f8').tobytes() for d in op['rotations'])
        residuals=[];nominal_residuals=[]
        for i,a in enumerate(structure.atoms):
            transformed=t.operation.rotation_frac@a.frac+t.operation.translation_frac
            target=structure.atoms[op['perm'][i]].frac+op['image_shift'][i]
            residuals.append(float(np.linalg.norm(cell.T@(transformed-target))))
            nominal_target=cell.T@target
            nominal_source=op['rotation_cart_nominal']@a.cart+cell.T@t.operation.translation_frac
            nominal_residuals.append(float(np.linalg.norm(nominal_source-nominal_target)))
        operation_payload.append(dict(index=int(t.operation.index),W=op['W'].tolist(),perm=op['perm'].tolist(),
            image_shift=op['image_shift'].tolist(),rotation_cart=op['rotation_cart_nominal'].tolist(),
            rotation_cart_nominal=op['rotation_cart_nominal'].tolist(),rotation_cart_cell_derived=t.operation.rotation_cart.tolist(),
            nominal_operation_index=op['nominal_index'],
            nominal_vs_cell_cartesian_max_difference=float(np.max(abs(op['rotation_cart_nominal']-t.operation.rotation_cart))),
            nominal_vs_cell_orbital_max_difference=max(float(np.max(abs(a-b))) for a,b in zip(op['rotations'],op['cell_derived_rotations'])),
            cell_derived_orbital_rotation_species=[d.tolist() for d in op['cell_derived_rotations']],
            translation_frac=t.operation.translation_frac.tolist(),determinant=float(np.linalg.det(t.operation.rotation_cart)),
            max_geometry_residual_Ang=max(residuals),max_nominal_geometry_residual_Ang=max(nominal_residuals),
            max_orbital_orthogonality_error=max(float(np.max(abs(d.T@d-np.eye(len(d))))) for d in op['rotations']),
            max_cell_derived_orbital_orthogonality_error=max(float(np.max(abs(d.T@d-np.eye(len(d))))) for d in op['cell_derived_rotations']),
            orbital_rotation_species=[d.tolist() for d in op['rotations']]))
    output_path.parent.mkdir(parents=True,exist_ok=True)
    output_path.write_bytes(b''.join(chunks))
    relevant=[ANGULAR_MODULE,PACKAGE_ROOT/'comp2real.json',Path(__file__).resolve()]
    approximate=[dict(g,rotations=g['cell_derived_rotations']) for g in operations]
    approximate_error=group_composition_error(approximate,orbital_tol=1e-6)
    geometry_hash=hashlib.sha256(cell.astype('<f8').tobytes()+xyz.astype('<f8').tobytes()+atom_species.astype('<i4').tobytes()+species_norb.astype('<i4').tobytes()).hexdigest()
    provenance=dict(schema='RSD3 little-endian v1; row-major; no padding',operation_count=len(operations),input=str(input_path),input_sha256=sha256(input_path),
        output=str(output_path),output_sha256=sha256(output_path),output_bytes=output_path.stat().st_size,
        natoms=len(xyz),N=structure.nwann,nspecies=len(species),species=species,species_basis=structure.species_basis,
        species_norb=species_norb.tolist(),atom_species=atom_species.tolist(),geometry_sha256=geometry_hash,
        geometry_modified=False,image_composition='relative shifts close exactly; common lattice translations cancel from every neighbor action',representation='real localized OpenMX NSOC orbitals; no spin factors; generalized H/S',
        package_root=str(PACKAGE_ROOT),package_module=str(ANGULAR_MODULE),
        package_hashes={str(p):sha256(p) for p in relevant},generator=str(Path(__file__).resolve()),generator_sha256=sha256(__file__),
        nominal_geometry_helper=str(NOMINAL_HELPER),nominal_geometry_helper_sha256=sha256(NOMINAL_HELPER),
        symprec_Ang=symprec,atom_tol_fractional=atom_tol,max_group_composition_error=err,operations=operation_payload,
        orbital_acceptance_bound=1e-12,max_cell_derived_group_composition_error=approximate_error,
        orbital_representation_exactness=('Exact nominal D3 Q3/Q2 from the established geometry helper for the six-operation D3 case; '
          'for other scalar groups the repository-local validated OpenMX AO blocks are used. Validated W/permutation/image mappings '
          'and atomic geometry are unchanged. Cell-derived package angular differences are separately recorded.'),
        environment_prefix=sys.prefix,python_executable=sys.executable,python_version=platform.python_version(),
        hostname=socket.gethostname(),working_directory=os.getcwd(),command_argv=sys.argv,
        process_threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')},
        category='metadata and mathematical validation; no native SCF result',final_process_state='metadata generation succeeded')
    provenance_path=Path(provenance_path) if provenance_path else output_path.with_suffix('.json')
    provenance_path.write_text(json.dumps(provenance,indent=2)+'\n')
    provenance['wall_seconds']=time.perf_counter()-started
    provenance['wall_timer_scope']='generate entry through binary and initial provenance JSON write; excludes only final timing annotation rewrite and process imports/startup'
    provenance['module_import_seconds']=MODULE_IMPORT_SECONDS
    provenance['module_import_timer_scope']='numpy and pinned package imports; shared once if generate called repeatedly in one process'
    provenance_path.write_text(json.dumps(provenance,indent=2)+'\n')
    return provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--provenance',type=Path);p.add_argument('--symprec',type=float,default=1e-4)
    p.add_argument('--atom-tol',type=float,default=5e-6)
    p.add_argument('--all-operations',action='store_true',help='include improper operations such as mirrors and inversion')
    a=p.parse_args(); result=generate(a.input,a.output,provenance_path=a.provenance,symprec=a.symprec,atom_tol=a.atom_tol,include_improper=a.all_operations)
    print(json.dumps({k:result[k] for k in ('input','output','natoms','N','species_norb','output_bytes','output_sha256','max_group_composition_error')},indent=2))

if __name__=='__main__':main()
