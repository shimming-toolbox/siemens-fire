import cupy as cp
import math

class hankel:
    """Basic Hankel Matrix Class"""

    @staticmethod
    def fwd(x, kernel):

        # Reshape to have coil dimension as dim2
        #x = np.reshape(x, (x.shape[0], x.shape[1], -1))
        nc = x.shape[2]
        
        # Get kernel-reduced dimensions
        #dimr = tuple(map(lambda i,j:i-j+1, x.shape[:2], kernel))
        dimr = (x.shape[0]-kernel[0]+1, x.shape[1]-kernel[1]+1)
        
        # Initialise matrix components
        h = cp.zeros((math.prod(dimr), math.prod(kernel), nc), dtype=cp.complex64)
        
        # Loop over all kernel locs
        for kx in range(dimr[0]):
            for ky in range(dimr[1]):
                h[kx*dimr[1]+ky,...] = x[kx:kx+kernel[0],ky:ky+kernel[1],:].copy().reshape((1,-1,nc))
                
        return h

    @staticmethod
    def adj(h, dims, kernel):

        # get coil dimension
        nc = h.shape[2]
        
        # get kernel-reduced dimensions
        #dimr = tuple(map(lambda i,j:i-j+1, dims, kernel))
        dimr = (dims[0]-kernel[0]+1, dims[1]-kernel[1]+1)

        # initialise output
        x = cp.zeros((dims[0], dims[1], nc), dtype=cp.complex64)
        
        # loop over all kernel locs
        for kx in range(dimr[0]):
            for ky in range(dimr[1]):
                x[kx:kx+kernel[0],ky:ky+kernel[1],:] += h[kx*dimr[1]+ky,...].reshape((kernel[0], kernel[1], nc))
        
        return x
    
    @staticmethod
    def norm(dims, kernel):
        
        if dims[0] >  2*kernel[0]:
            U = cp.concatenate((cp.arange(kernel[0])+1, kernel[0]*cp.ones((dims[0]-2*kernel[0])), cp.arange(kernel[0],0,-1)))
        elif dims[0] > kernel[0]:
            U = cp.arange((dims[0]+1)//2)+1
            U = cp.concatenate((U[:(dims[0]+1)//2],cp.flip(U[:(dims[0])//2])))
        else:
            U = cp.ones((kernel[0],))
        if dims[1] > 2*kernel[1]:
            V = cp.concatenate((cp.arange(kernel[1])+1, kernel[1]*cp.ones((dims[1]-2*kernel[1])), cp.arange(kernel[1],0,-1)))
        elif dims[1] > kernel[1]:
            V = cp.arange((dims[1]+1)//2)+1
            V = cp.concatenate((V[:(dims[0]+1)//2],cp.flip(V[:(dims[0])//2])))
        else:
            V = cp.ones((kernel[1],))
            
        N = U[:,None]@V[None,:]
        
        return N[:,:,None]

    @staticmethod
    def size(dims, kernel):
        return ((dims[0]-kernel[0]+1)*(dims[1]-kernel[1]+1), math.prod(kernel))
    
class c_matrix:
    """LORAKS-style C-Matrix"""
    
    @staticmethod
    def fwd(x, kernel):

        c = hankel.fwd(x, kernel)
        
        return c.reshape((c.shape[0],-1))
        
    @staticmethod
    def adj(c, dims, kernel):
        
        c = c.reshape(c_matrix.size(dims,kernel)+(-1,))
        
        return hankel.adj(c, dims, kernel)

    @staticmethod
    def norm(dims, kernel):
        return hankel.norm(dims, kernel)

    @staticmethod
    def size(dims, kernel):
        return hankel.size(dims, kernel)
    

class s_matrix:
    """LORAKS-style S-Matrix"""
    
    @staticmethod
    def fwd(x, kernel):

        s_pos = hankel.fwd(x, kernel)
        s_neg = cp.flip(s_pos, axis=0)
        
        s = cp.concatenate((cp.concatenate((cp.real(s_pos)-cp.real(s_neg), 
                                            cp.imag(s_neg)-cp.imag(s_pos)),axis=1),
                            cp.concatenate((cp.imag(s_pos)+cp.imag(s_neg), 
                                            cp.real(s_pos)+cp.real(s_neg)),axis=1)),axis=0)
        return s.reshape((s.shape[0],-1))

        
    @staticmethod
    def adj(s, dims, kernel):

        s = s.reshape(s_matrix.size(dims,kernel)+(-1,))

        A = s[:s.shape[0]//2, :s.shape[1]//2, :]
        B = s[:s.shape[0]//2, s.shape[1]//2:, :]
        C = s[s.shape[0]//2:, :s.shape[1]//2, :]
        D = s[s.shape[0]//2:, s.shape[1]//2:, :]

        s_pos = A+D + 1j*(C-B)
        s_neg = D-A + 1j*(B+C)

        s_pos = s_pos + cp.flip(s_neg, axis=0)

        return hankel.adj(s_pos, dims, kernel)
    
    @staticmethod
    def norm(dims, kernel):
        return 4*c_matrix.norm(dims, kernel)
    
    @staticmethod
    def size(dims, kernel):
        return c_matrix.size(dims, kernel)*cp.array([2,2])

class vcc_matrix:
    """Virtual conjugate channel matrix, alternative formulation of S-Matrix"""
    
    @staticmethod
    def fwd(x, kernel):

        v_pos = hankel.fwd(x, kernel)
        v_neg = cp.conj(cp.flip(v_pos))
        
        v = cp.concatenate((v_pos, v_neg), axis=1)
        
        return v.reshape((v.shape[0],-1))
        
    @staticmethod
    def adj(v, dims, kernel):
        
        v = v.reshape(vcc_matrix.size(dims,kernel)+(-1,))

        v_pos = v[:, :v.shape[1]//2, :]
        v_neg = v[:, v.shape[1]//2:, :]

        v_pos = v_pos + cp.flip(cp.conj(v_neg))
        
        return hankel.adj(v_pos, dims, kernel)
    
    @staticmethod
    def norm(dims, kernel):
        return 2*hankel.norm(dims, kernel)
    
    @staticmethod
    def size(dims, kernel):
        return tuple(x * y for x, y in zip(hankel.size(dims, kernel), (1,2)))
    
        
def ADMM(d, mtx, kernel, r, p=1E-2, niters=100, tol=1E-4, init=None, mode='hard'):
    """ADMM reconstruction for Hankel-structured low-rank optimization.
       If mode == 'hard', parameter r is interpreted as an integer rank constraint,
       and this solves the non-convex strict rank constraint optimization problem
       If mode == 'soft', parameter r is interpreted as a scalar lambda weighting 
       on the nuclear norm penalty, and solves the convex nuclear norm minimization"""

    # Check threshold type
    if mode == 'hard':
        r = int(r)
    elif mode == 'soft':
        lam = float(r) 
    else:
        raise ValueError(f'Unsupported mode: {mode!r}')
	
    
    # Pad input
    d, crop = __pad(cp.asarray(d, dtype=cp.complex64))
    dims = d.shape[:2]
    
    # Get sampling mask
    M = (d !=0)
    
    # Initialise 
    if init is None:
        init = d
    else:
        init, _ = __pad(cp.asarray(init, dtype=cp.complex64))
        
    x = init
    z = mtx.fwd(x, kernel)
    u = 0*z
    
    # Get normalisation factor
    N = mtx.norm(dims, kernel)
    
    # Precompute LHS
    Q = 1/(M + (p/2)*N)
    Q[cp.isinf(Q)] = 0
    
    # ADMM iterations
    for i in range(niters):

        print(f'\rIter {i:04d}', end='\r')
        
        # x-update
        xx = Q*(d + (p/2)*mtx.adj(z - u, dims, kernel))

        # z-update
        H = mtx.fwd(xx, kernel)

        if mode == 'hard': 
                S,V = half_SVD(H + u)
                z = ((H + u)@V[:,:r])@cp.conj(V[:,:r].T)
                last_sv = float(S[r - 1]) if r <= len(S) else float(S[-1])
        elif mode == 'soft':
                S,V = half_SVD(H + u)
                mask = S > lam
                S2 = S[mask] - lam
                z = (((H + u)@V[:,mask])*(S2/S[mask]))@ cp.conj(V[:,mask].T)
                last_sv = float(S[mask][-1]) if cp.sum(mask) > 0 else 0.0
        
        # u-update
        u = u + H - z
        
        # Check relative update tolerance
        update = cp.linalg.norm(xx.ravel()-x.ravel())/cp.linalg.norm(x.ravel())
        if update < tol and i > 0:
            print(f'Min Update Tolerance Reached at {i} iterations')
            break
        
        # Save estimate
        x = xx 
    print(f'Last singular value kept at exit: {last_sv:.6e}')

    return cp.asnumpy(__unpad(x, crop))
    

    
## Helper
def __pad(d):
    # Zero-pad input if even, to make k-space symmetric about origin
    # Also add coil dimension if not present
    if cp.ndim(d) == 2:
        d = d[:,:,None]
        
    #pad = 1 - cp.mod(d.shape[:2],2)
    pad = (1-d.shape[0]%2, 1-d.shape[1]%2)
    d_cpu = cp.asnumpy(d)
    d_padded_cpu = cp.pad(d_cpu, ((0, pad[0]), (0, pad[1]), (0, 0)), mode='wrap')
    d = cp.asarray(d_padded_cpu)
    crop = (d.shape[0] - pad[0], d.shape[1] - pad[1])
    
    return d, crop

def __unpad(x, crop):
    if cp.ndim(x) == 2:
        x = x[:,:,None]
        
    if x.shape[2] == 1:
        return x[:crop[0],:crop[1],0]
    else:
        return x[:crop[0],:crop[1],:]


def half_SVD(x):

    d, v = cp.linalg.eigh(cp.conj(x).T@x)
    
    ii = cp.argsort(cp.abs(d))
    
    s = cp.sqrt(d[ii[::-1]])
    v = v[:,ii[::-1]] 
    
    return s, v    
