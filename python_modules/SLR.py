import numpy as np
from matplotlib import pyplot as plt
import utils
'''
This script was made by Prof. Mark Chiew at the University of Toronto.
'''
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
        h = np.zeros((np.prod(np.array(dimr)), np.prod(np.array(kernel)), nc), dtype=np.complex64)
        
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
        x = np.zeros((dims[0], dims[1], nc), dtype=np.complex64)
        
        # loop over all kernel locs
        for kx in range(dimr[0]):
            for ky in range(dimr[1]):
                x[kx:kx+kernel[0],ky:ky+kernel[1],:] += h[kx*dimr[1]+ky,...].reshape((kernel[0], kernel[1], nc))
        
        return x
    
    @staticmethod
    def norm(dims, kernel):
        
        if dims[0] >  2*kernel[0]:
            U = np.concatenate((np.arange(kernel[0])+1, kernel[0]*np.ones((dims[0]-2*kernel[0])), np.arange(kernel[0],0,-1)))
        elif dims[0] > kernel[0]:
            U = np.arange((dims[0]+1)//2)+1
            U = np.concatenate((U[:(dims[0]+1)//2],np.flip(U[:(dims[0])//2])))
        else:
            U = np.ones((kernel[0],))
        if dims[1] > 2*kernel[1]:
            V = np.concatenate((np.arange(kernel[1])+1, kernel[1]*np.ones((dims[1]-2*kernel[1])), np.arange(kernel[1],0,-1)))
        elif dims[1] > kernel[1]:
            V = np.arange((dims[1]+1)//2)+1
            V = np.concatenate((V[:(dims[0]+1)//2],np.flip(V[:(dims[0])//2])))
        else:
            V = np.ones((kernel[1],))
            
        N = U[:,np.newaxis]@V[np.newaxis,:]
        
        return N[:,:,np.newaxis]

    @staticmethod
    def size(dims, kernel):
        return ((dims[0]-kernel[0]+1)*(dims[1]-kernel[1]+1), np.prod(kernel))
    
class c_matrix:
    """LORAKS-style C-Matrix"""
    
    @staticmethod
    def fwd(x, kernel):

        c = hankel.fwd(x, kernel)
        
        return c.reshape((c.shape[0],-1))
        
    @staticmethod
    def adj(c, dims, kernel):
        
        c = c.reshape(np.concatenate((c_matrix.size(dims,kernel),(-1,))))
        
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
        s_neg = np.flip(s_pos, axis=0)
        
        s = np.concatenate((np.concatenate((np.real(s_pos)-np.real(s_neg), 
                                            np.imag(s_neg)-np.imag(s_pos)),axis=1),
                            np.concatenate((np.imag(s_pos)+np.imag(s_neg), 
                                            np.real(s_pos)+np.real(s_neg)),axis=1)),axis=0)
        return s.reshape((s.shape[0],-1))

        
    @staticmethod
    def adj(s, dims, kernel):

        s = s.reshape(np.concatenate((s_matrix.size(dims,kernel),(-1,))))

        A = s[:s.shape[0]//2, :s.shape[1]//2, :]
        B = s[:s.shape[0]//2, s.shape[1]//2:, :]
        C = s[s.shape[0]//2:, :s.shape[1]//2, :]
        D = s[s.shape[0]//2:, s.shape[1]//2:, :]

        s_pos = A+D + 1j*(C-B)
        s_neg = D-A + 1j*(B+C)

        s_pos = s_pos + np.flip(s_neg, axis=0)

        return hankel.adj(s_pos, dims, kernel)
    
    @staticmethod
    def norm(dims, kernel):
        return 4*c_matrix.norm(dims, kernel)
    
    @staticmethod
    def size(dims, kernel):
        return c_matrix.size(dims, kernel)*np.array([2,2])

class vcc_matrix:
    """Virtual conjugate channel matrix, alternative formulation of S-Matrix"""
    
    @staticmethod
    def fwd(x, kernel):

        v_pos = hankel.fwd(x, kernel)
        v_neg = np.conj(np.flip(v_pos))
        
        v = np.concatenate((v_pos, v_neg), axis=1)
        
        return v.reshape((v.shape[0],-1))
        
    @staticmethod
    def adj(v, dims, kernel):
        
        v = v.reshape(np.concatenate((vcc_matrix.size(dims,kernel),(-1,))))

        v_pos = v[:, :v.shape[1]//2, :]
        v_neg = v[:, v.shape[1]//2:, :]

        v_pos = v_pos + np.flip(np.conj(v_neg))
        
        return hankel.adj(v_pos, dims, kernel)
    
    @staticmethod
    def norm(dims, kernel):
        return 2*hankel.norm(dims, kernel)
    
    @staticmethod
    def size(dims, kernel):
        return hankel.size(dims, kernel)*np.array([1,2])
    
        
def ADMM(d, mtx, kernel, r, p=1E-2, niters=100, tol=1E-4, init=None):
    """Non-convex ADMM reconstruction with strict rank constraints"""
    
    # Pad input
    d, crop = __pad(d)
    dims = d.shape[:2]
    
    # Get sampling mask
    M = (d !=0)
    
    # Initialise 
    if init is None:
        init = d
    else:
        init, _ = __pad(init)
        
    x = init
    z = mtx.fwd(x, kernel)
    u = 0*z
    
    # Get normalisation factor
    N = mtx.norm(dims, kernel)
    
    # Precompute LHS
    Q = 1/(M + (p/2)*N)
    Q[np.isinf(Q)] = 0
    
    # ADMM iterations
    for i in range(niters):

        print(f'\rIter {i:04d}', end='\r')
        
        # x-update
        xx = Q*(d + (p/2)*mtx.adj(z - u, dims, kernel))
        
        # z-update
        H = mtx.fwd(xx, kernel)
        _,V = half_SVD(H + u)
        z = (H + u)@(V[:,:r]@np.conj(V[:,:r].T))
        
        # u-update
        u = u + H - z
        
        # Check relative update tolerance
        update = np.linalg.norm(xx.ravel()-x.ravel())/np.linalg.norm(x.ravel())
        if update < tol and i > 0:
            print(f'Min Update Tolerance Reached at {i} iterations')
            break
        
        # Save estimate
        x = xx 

    return __unpad(x, crop)
    
    
## Helpers
def __pad(d):
    # Zero-pad input if even, to make k-space symmetric about origin
    # Also add coil dimension if not present
    if np.ndim(d) == 2:
        d = d[:,:,np.newaxis]
        
    pad = 1 - np.mod(d.shape[:2],2)
    d = np.pad(d, ((0,pad[0]),(0,pad[1]),(0,0)), mode='wrap')
    crop = d.shape[:2] - pad
    
    return d, crop

def __unpad(x, crop):
    if np.ndim(x) == 2:
        x = x[:,:,np.newaxis]
        
    if x.shape[2] == 1:
        return x[:crop[0],:crop[1],0]
    else:
        return x[:crop[0],:crop[1],:]

def half_SVD(x):

    d, v = np.linalg.eigh(np.conj(x).T@x)
    
    ii = np.argsort(np.abs(d))
    
    s = np.sqrt(d[ii[::-1]])
    v = v[:,ii[::-1]] 
    
    return s, v
