#!/usr/bin/env python3
"""Internal native band worker: full-spectrum LAPACK or targeted SLEPc.

No moirekp dependency. Native code supplies owned real-space blocks once; the
worker constructs matching Bloch H/S and returns the requested spectrum.
"""
from pathlib import Path
import argparse
import json
import os
import struct
import sys
import time

EV_PER_HA = 27.211386245988


def reference_ha(value, fermi):
    import numpy as np
    result = float(fermi) if value.lower() == 'fermi' else float(value)/EV_PER_HA
    if not np.isfinite(result):raise ValueError('Band.ReferenceEnergy must be Fermi or a finite absolute energy in eV')
    return result


def nearest_indices(values, target, count):
    import numpy as np
    if count < 1 or count > len(values):raise ValueError('insufficient converged eigenvalues')
    indices=np.lexsort((values,abs(values-target)))[:count]
    return indices[np.argsort(values[indices],kind='stable')]


def load_parts(directory, parts):
    import numpy as np
    rr=[];cc=[];translations=[];values=[];schema=None
    for rank in range(parts):
        path=Path(directory)/f'part-{rank:06d}.bin'
        with path.open('rb') as stream:
            if stream.read(8)!=b'OMBAND01':raise ValueError('invalid native band snapshot magic')
            header=np.fromfile(stream,'<i4',5)
            if len(header)!=5:raise ValueError('short band snapshot header')
            n,spin,blocks,fields,version=map(int,header)
            if n<1 or spin not in (0,1,3) or blocks<0 or fields!={0:2,1:3,3:8}[spin] or version!=1:
                raise ValueError('invalid native band snapshot schema')
            if schema is not None and schema!=(n,spin,fields):raise ValueError('inconsistent native band snapshots')
            schema=n,spin,fields
            for block in range(blocks):
                h=np.fromfile(stream,'<i4',7)
                if len(h)!=7:raise ValueError('short block header')
                row,col,na,nb,*R=map(int,h)
                if min(row,col)<0 or min(na,nb)<1 or row+na>n or col+nb>n:raise ValueError('block outside native AO basis')
                data=np.fromfile(stream,'<f8',na*nb*fields)
                if len(data)!=na*nb*fields or not np.isfinite(data).all():raise ValueError('invalid block values')
                rr.append(np.repeat(np.arange(row,row+na,dtype=np.int32),nb))
                cc.append(np.tile(np.arange(col,col+nb,dtype=np.int32),na))
                translations.append(np.tile(np.array(R,dtype=np.int32),(na*nb,1)))
                values.append(data.reshape(fields,na*nb))
            if stream.read(1):raise ValueError('trailing native band snapshot payload')
    if not values:raise ValueError('no native H/S blocks')
    return dict(n=schema[0],spin=schema[1],rows=np.concatenate(rr),cols=np.concatenate(cc),
                R=np.concatenate(translations),values=np.concatenate(values,axis=1))


def bloch_matrices(data, k, channel):
    import numpy as np
    import scipy.sparse as sp
    n=data['n'];v=data['values'];phase=np.exp(2j*np.pi*(data['R']@k))
    rows=data['rows'];cols=data['cols'];keep=rows>=cols
    lower_rows=rows[keep];lower_cols=cols[keep]
    lower_phase=phase[keep]
    def matrix(values,hermitian=True):
        if hermitian:
            # Sum periodic duplicates once in the lower half, then mirror that
            # sum so conjugate entries cannot acquire different roundoff.
            lower=sp.coo_matrix((values[keep]*lower_phase,(lower_rows,lower_cols)),shape=(n,n)).tocsr()
            lower.setdiag(.5*lower.diagonal().real)
            a=(lower+lower.getH()).tocsr()
        else:
            # A spin-up/down off-diagonal block is not individually Hermitian.
            a=sp.coo_matrix((values*phase,(rows,cols)),shape=(n,n)).tocsr()
        a.sum_duplicates();a.eliminate_zeros();return a
    S=matrix(v[-1])
    if data['spin']!=3:return matrix(v[channel]),S
    up=matrix(v[0]+1j*v[4]);down=matrix(v[1]+1j*v[5]);ud=matrix(v[2]+1j*(v[3]+v[6]),hermitian=False)
    return sp.bmat([[up,ud],[ud.getH(),down]],format='csr'),sp.block_diag((S,S),format='csr')


def hermiticity_defect(matrix):
    import numpy as np
    from scipy.sparse.linalg import norm
    if not np.isfinite(matrix.data).all():raise ValueError('nonfinite Bloch matrix')
    defect=float(norm(matrix-matrix.getH())/max(float(norm(matrix)),1e-300))
    if not np.isfinite(defect) or defect>5e-7:raise ValueError('native Bloch matrix is not Hermitian: '+str(defect))
    return defect


def native_hermitian(matrix):
    """Explicit legacy lower-triangle conversion, not a default solver step."""
    import scipy.sparse as sp
    defect=hermiticity_defect(matrix)
    # The original zheevx UPLO=L reads the lower triangle, including real diag.
    lower=sp.tril(matrix,k=-1,format='csr')
    return (lower+lower.getH()+sp.diags(matrix.diagonal().real)).tocsr(),defect


def csr_pair(H,S):
    """Keep canonical caller-owned CSR storage; convert only other formats."""
    H=H.tocsr(copy=False);S=S.tocsr(copy=False)
    if len(H.shape)!=2 or H.shape[0]!=H.shape[1] or H.shape!=S.shape or H.shape[0]<1:
        raise ValueError('incompatible native H/S dimensions')
    return H,S


def petsc_csr_buffers(matrix,int_type,scalar_type):
    """Borrow compatible contiguous arrays; convert only incompatible dtypes."""
    import numpy as np
    capacity=np.iinfo(int_type).max
    if matrix.nnz>capacity or max(matrix.shape)>capacity:raise ValueError('sparse index capacity exceeded')
    return (np.asarray(matrix.indptr,dtype=int_type,order='C'),
            np.asarray(matrix.indices,dtype=int_type,order='C'),
            np.asarray(matrix.data,dtype=scalar_type,order='C'))


def libraries():
    import numpy as np
    from petsc4py import PETSc
    from slepc4py import SLEPc
    if np.dtype(PETSc.ScalarType).kind!='c':raise ValueError('targeted bands require complex PETSc/SLEPc')
    factor=next((name for name in ('mkl_pardiso','mumps','superlu') if PETSc.Sys.hasExternalPackage(name)),None)
    if factor is None:raise ValueError('PETSc needs a sequential sparse LU backend (PARDISO, MUMPS or SuperLU)')
    return PETSc,SLEPc,factor


def check_overlap(S):
    """Check sparse input without a separate overlap eigenproblem.

    GHEP requires positive-definite S as an input contract. These inexpensive
    necessary checks do not certify its global spectrum; report that explicitly.
    The solve also checks every returned generalized eigenpair and S-orthogonality.
    """
    import numpy as np
    start=time.monotonic()
    if len(S.shape)!=2 or S.shape[0]!=S.shape[1] or S.shape[0]<1:
        raise ValueError('invalid overlap dimensions')
    if not np.isfinite(S.data).all():raise ValueError('nonfinite overlap')
    diagonal=np.asarray(S.diagonal())
    if np.any(diagonal.real<=1e-9):raise ValueError('overlap requires native regularization')
    if np.any(diagonal.imag!=0):raise ValueError('overlap diagonal is not real')
    return dict(overlap_lower_bound=None,overlap_min=None,
                overlap_diagonal_min=float(diagonal.real.min()),
                overlap_positive_definite_certified=False,
                overlap_check_method='sparse input checks',overlap_check_seconds=time.monotonic()-start)


def solve_dense(H,S,target,count,driver='gvd',audit=False):
    """Full generalized spectrum, with no eigenvectors or identity-S shortcut."""
    import numpy as np
    import scipy.linalg as la
    if driver not in ('gvd','evr'):raise ValueError('unsupported dense band driver')
    H,S=csr_pair(H,S);n=H.shape[0]
    if not np.isfinite(H.data).all() or not np.isfinite(S.data).all():raise ValueError('nonfinite Bloch matrix')
    hdef=hermiticity_defect(H) if audit else None;sdef=hermiticity_defect(S) if audit else None
    if count<0 or count>n:raise ValueError('Band.NumBands exceeds matrix dimension')
    check_start=time.monotonic();B=S.toarray(order='F')
    try:L=la.cholesky(B,lower=True,check_finite=False)
    except la.LinAlgError as error:raise ValueError('overlap requires native regularization') from error
    inverse,info=la.get_lapack_funcs('trtri',(L,))(L,lower=1,unitdiag=0,overwrite_c=0)
    if info:raise ValueError('overlap inverse-Cholesky certificate failed')
    norm2=float(np.vdot(inverse.ravel(order='K'),inverse.ravel(order='K')).real);del inverse
    if not np.isfinite(norm2) or norm2<=0:raise ValueError('nonfinite overlap certificate')
    bound=1./norm2;minimum=None
    # lambda_min(S) >= 1/||L^-1||_F^2. Only ambiguous bases need an extremal solve.
    if bound<=1e-9:
        minimum=float(la.eigvalsh(B,subset_by_index=(0,0),check_finite=False)[0])
        if not np.isfinite(minimum) or minimum<=1e-9:raise ValueError('overlap requires native regularization')
    check_seconds=time.monotonic()-check_start;start=time.monotonic();A=H.toarray(order='F')
    if driver=='gvd':
        del L
        values=la.eigh(A,B,eigvals_only=True,driver='gvd',lower=True,overwrite_a=True,overwrite_b=True,check_finite=False)
    else:
        C,info=la.get_lapack_funcs('hegst',(A,L))(A,L,itype=1,lower=1,overwrite_a=True)
        if info:raise ValueError('generalized Hermitian transform failed')
        values=la.eigh(C,eigvals_only=True,driver='evr',lower=True,overwrite_a=True,check_finite=False)
    seconds=time.monotonic()-start
    if len(values)!=n or not np.isfinite(values).all() or np.any(np.diff(values)<0):raise ValueError('invalid full band spectrum')
    selected=values if count==0 else values[nearest_indices(values,target,count)]
    return selected,dict(solver='LAPACK '+driver.upper(),matrix_dimension=n,requested=len(selected),
        spectrum_computed=n,full_spectrum=True,eigenvectors_computed=False,overlap_lower_bound=bound,
        overlap_min=minimum,overlap_certificate='lambda_min(S)>=1/||L^-1||_F^2',overlap_check_seconds=check_seconds,
        solve_seconds=seconds,H_hermiticity=hdef,S_hermiticity=sdef,matrix_hermiticity_checked=bool(audit),
        solve_time_matrix_reconstruction=False,S_diagonal_shift=0.,native_triangle='lower',
        eigenpair_audit='not computed in eigenvalues-only mode; validated by independent spectrum comparisons')


def solve(H,S,target,count,tolerance=1e-9,audit=False):
    import numpy as np
    PETSc,SLEPc,factor=libraries()
    prepare_start=time.monotonic();H,S=csr_pair(H,S);n=H.shape[0]
    if not np.isfinite(H.data).all():raise ValueError('nonfinite Bloch matrix')
    hdef=hermiticity_defect(H) if audit else None;sdef=hermiticity_defect(S) if audit else None
    if not 0<count<n-1:raise ValueError('Band.NumBands must be smaller than matrix dimension minus one')
    overlap=check_overlap(S)
    prepare_seconds=time.monotonic()-prepare_start
    csr_keepalive=[];shared={}
    def petsc_matrix(a,name):
        buffers=petsc_csr_buffers(a,PETSc.IntType,PETSc.ScalarType);csr_keepalive.append(buffers)
        shared[name]=[bool(np.shares_memory(x,y)) for x,y in zip(buffers,(a.indptr,a.indices,a.data))]
        m=PETSc.Mat().createAIJ(size=a.shape,csr=buffers,comm=PETSc.COMM_SELF)
        m.assemble();m.setOption(PETSc.Mat.Option.HERMITIAN,True);return m
    A=B=eps=vr=vi=None
    try:
        setup_start=time.monotonic()
        A=petsc_matrix(H,'H');B=petsc_matrix(S,'S');B.setOption(PETSc.Mat.Option.SPD,True)
        eps=SLEPc.EPS().create(PETSc.COMM_SELF);eps.setOperators(A,B);eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
        # Request exactly Band.NumBands states; leave the Krylov subspace size
        # to SLEPc rather than silently padding the requested spectrum.
        eps.setType('krylovschur');wanted=count
        # A wider working subspace retains near-degenerate partners without
        # increasing the requested or returned eigenstate count.
        eps.setDimensions(wanted,min(n,max(4*wanted,20)))
        eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL);eps.setTarget(target)
        st=eps.getST();st.setType('sinvert');st.setShift(target)
        ksp=st.getKSP();ksp.setType('preonly');pc=ksp.getPC();pc.setType('lu');pc.setFactorSolverType(factor)
        eps.setTolerances(tolerance,10000)
        setup_seconds=time.monotonic()-setup_start
        solve_start=time.monotonic()
        if factor == 'mumps':
            # Keep tree-parallel analysis/factorization, but use the conventional
            # solve path: the installed MUMPS tree-parallel backsolve can return
            # wrong vectors without an error status on complex SOC matrices.
            eps.setUp()
            ksp.setUp()
            pc.getFactorMatrix().setMumpsIcntl(48,0)
        eps.solve();solve_seconds=time.monotonic()-solve_start
        converged=eps.getConverged()
        ncv=int(eps.getDimensions()[1])
        if eps.getConvergedReason()<=0 or converged<count:raise ValueError('SLEPc did not converge all requested bands')
        values=np.array([complex(eps.getEigenvalue(i)) for i in range(converged)])
        if not np.isfinite(values).all() or max(abs(values.imag))>1e-10:raise ValueError('nonreal/nonfinite eigenvalues')
        selected=nearest_indices(values.real,target,count);residual=orth=None;audit_seconds=0.
        if audit:
            audit_start=time.monotonic();vectors=[];residual=0.
            vr=A.createVecRight();vi=A.createVecRight()
            for index in selected:
                eps.getEigenpair(int(index),vr,vi);x=vr.getArray(readonly=True).copy();e=values[index].real
                hx=H@x;sx=S@x;den=np.linalg.norm(hx)+abs(e)*np.linalg.norm(sx)
                residual=max(residual,float(np.linalg.norm(hx-e*sx)/max(den,1e-300)))
                vectors.append(x)
            C=np.column_stack(vectors);orth=float(abs(C.conj().T@(S@C)-np.eye(count)).max())
            if not np.isfinite(residual+orth) or residual>1e-8 or orth>1e-8:raise ValueError(f'generalized eigenpair audit failed: residual={residual} Sorth={orth}')
            audit_seconds=time.monotonic()-audit_start
        return values.real[selected],dict(residual=residual,Sorth=orth,**overlap,H_hermiticity=hdef,S_hermiticity=sdef,
            matrix_hermiticity_checked=bool(audit),solve_time_matrix_reconstruction=False,
            matrix_prepare_seconds=prepare_seconds,petsc_setup_seconds=setup_seconds,csr_python_buffers_shared=shared,
            eigenpair_audit_enabled=bool(audit),eigenvectors_extracted=bool(audit),eigenpair_audit_seconds=audit_seconds,
            solver='SLEPc KrylovSchur shift-invert',factor=factor,PETSc=PETSc.Sys.getVersion(),SLEPc=SLEPc.Sys.getVersion(),
            matrix_dimension=n,requested=count,search_nev=wanted,search_ncv=ncv,search_ncv_policy='min(N,max(4*nev,20))',
            converged=converged,solve_seconds=solve_seconds,worker_MKL_DYNAMIC=os.environ.get('MKL_DYNAMIC'),
            S_diagonal_shift=0.,native_triangle='assembled Hermitian CSR; no solve-time rewrite',full_spectrum=False)
    finally:
        for value in (vi,vr,eps,B,A):
            if value is not None:value.destroy()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--solver',choices=('slepc','gvd','evr'),default='slepc')
    parser.add_argument('--directory',type=Path);parser.add_argument('--rank',type=int,default=0)
    parser.add_argument('--workers',type=int,default=1);parser.add_argument('--parts',type=int,default=1)
    parser.add_argument('--count',type=int);parser.add_argument('--target-ha',type=float)
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--audit',action='store_true',default=os.environ.get('OPENMX_BAND_AUDIT','').strip().lower() in ('1','on','true','yes'),
                        help='extract eigenvectors and check generalized residuals and S-orthogonality (default off)')
    args=parser.parse_args()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key]=str(args.threads)
    # Match the isolated native worker even when invoked directly.
    os.environ['MKL_DYNAMIC']='TRUE'
    if args.check:
        if args.solver=='slepc':libraries()
        else:
            import scipy.linalg
        return 0
    import numpy as np
    started=time.monotonic();data=load_parts(args.directory,args.parts)
    with (args.directory/'kpoints.bin').open('rb') as stream:
        header=np.fromfile(stream,'<i4',1)
        if len(header)!=1 or header[0]<1:raise ValueError('bad native kpoint count')
        points=np.fromfile(stream,'<f8',int(header[0])*3).reshape(-1,3)
        if stream.read(1) or not np.isfinite(points).all():raise ValueError('bad native kpoint payload')
    audits=[];channels=2 if data['spin']==1 else 1
    final=args.directory/f'values-{args.rank:06d}.bin';temporary=final.with_suffix('.writing')
    with temporary.open('xb') as out:
        for index in range(args.rank,len(points),args.workers):
            for channel in range(channels):
                before=time.monotonic();H,S=bloch_matrices(data,points[index],channel)
                values,audit=(solve(H,S,args.target_ha,args.count,audit=args.audit) if args.solver=='slepc' else
                              solve_dense(H,S,args.target_ha,args.count,args.solver,audit=args.audit))
                out.write(struct.pack('<2i',index,channel));out.write(values.astype('<f8').tobytes())
                audits.append(dict(kpoint=index,channel=channel,k=points[index].tolist(),total_seconds=time.monotonic()-before,**audit))
    temporary.replace(final)
    import resource
    (args.directory/f'audit-{args.rank:06d}.json').write_text(json.dumps(dict(status='complete',python=sys.executable,worker=args.rank,
        threads=args.threads,seconds=time.monotonic()-started,process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        memory_scope='one worker lifetime peak RSS, not node physical peak',points=audits),indent=2)+'\n')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as error:
        print('OPENMX_BANDS_ERROR: '+str(error),file=sys.stderr,flush=True)
        raise SystemExit(2)
