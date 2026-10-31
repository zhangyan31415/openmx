"""Exact scalar C3 conjugate-pair algebra; research/reference implementation."""
import numpy as np
from scipy.linalg import eigh


def realify(z):
    z = np.asarray(z)
    return np.block([[z.real, -z.imag], [z.imag, z.real]])


def paired_basis(rotation, tol=2e-11):
    """Return X,Y,A with R X=c X+s Y, R Y=-s X+c Y, R A=A.

    All matrices here are tiny atom-orbit/shell matrices, never full crystal U.
    The pivoted real construction fixes conjugate partners deterministically.
    """
    r = np.asarray(rotation, dtype=float)
    n = len(r)
    eye = np.eye(n)
    if r.shape != (n,n) or max(np.max(abs(r.T@r-eye), initial=0),
                              np.max(abs(r@r@r-eye), initial=0)) > tol:
        raise ValueError('rotation is not an orthogonal C3 representation')
    p0 = (eye+r+r.T)/3
    pe = eye-p0
    j = (r-r.T)/np.sqrt(3)
    if np.max(abs(j@j+pe), initial=0) > tol:
        raise ValueError('invalid C3 complex structure')
    # A separate projected coordinate pivot fixes the trivial subspace.
    def real_projector_basis(projector):
        work = projector.copy()
        out = []
        while True:
            norms = np.einsum('ij,ij->j', work, work)
            k = int(np.argmax(norms))
            if norms[k] < tol*tol:
                break
            v = work[:, k]/np.sqrt(norms[k])
            out.append(v)
            work -= np.outer(v, v@work)
        return np.column_stack(out) if out else np.empty((n,0))
    work = pe.copy()
    xs, ys = [], []
    while True:
        norms = np.einsum('ij,ij->j', work, work)
        k = int(np.argmax(norms))
        if norms[k] < tol*tol:
            break
        x = work[:,k]/np.sqrt(norms[k])
        y = j@x
        if abs(x@y) > tol or abs(y@y-1) > tol:
            raise ValueError('C3 conjugate partner is not orthonormal')
        xs.append(x)
        ys.append(y)
        work -= np.outer(x,x@work)+np.outer(y,y@work)
    x = np.column_stack(xs) if xs else np.empty((n,0))
    y = np.column_stack(ys) if ys else np.empty((n,0))
    a = real_projector_basis(p0)
    u = np.column_stack([x,y,a])
    if u.shape != (n,n) or np.max(abs(u.T@u-eye), initial=0) > tol:
        raise ValueError('paired basis does not span the local representation')
    return x,y,a


def _complexify(matrix, tol):
    n = len(matrix)
    if n%2:
        raise ValueError('C3 real E block must have even dimension')
    m = n//2
    xx, xy, yx, yy = matrix[:m,:m], matrix[:m,m:], matrix[m:,:m], matrix[m:,m:]
    # Average equivalent entries so tiny floating-point noise is not amplified.
    a, b = (xx+yy)/2, (yx-xy)/2
    z = a+1j*b
    scale = max(np.max(abs(matrix), initial=0), 1)
    if np.max(abs(matrix-realify(z)), initial=0)>tol*scale:
        raise ValueError('real block lacks the C3 complex structure')
    if np.max(abs(z-z.conj().T), initial=0)>tol*scale:
        raise ValueError('complex C3 block is not Hermitian')
    return z


def solve_pair(hr, sr, tol=2e-11):
    """Generalized m-complex solve, returning the full 2m-real eigenpairs."""
    hc, sc = _complexify(hr,tol), _complexify(sr,tol)
    e,z = eigh(hc,sc)
    m = len(e)
    c = np.empty((2*m,2*m))
    c[:,0::2] = np.vstack([z.real,z.imag])
    c[:,1::2] = np.vstack([-z.imag,z.real])
    return np.repeat(e,2),c
