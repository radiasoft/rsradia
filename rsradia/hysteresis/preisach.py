### preisach.py
### February 2023

from numpy import array, zeros, ceil, diff, sign, ndarray

from ..utils.distributions import *
from .hysteresis_model import HysteresisModel

class Preisach(HysteresisModel):
    '''
    A Preisach model of magnetic hysteresis

    Parameters:
        Ms: Saturation magnetization of the magnet
        ab_max: Maximum magnitude of hysteron on/off field values
        ab_res: Resolution of alpha/beta grid (point spacing)
        distribution: Type of weight distribution used
        dist_params: Parameters of the weight distribution
    '''
    
    # List slotted members of Preisach (non-extendable)
    __slots__ = ("_Ms", "_ab_max", "_ab_grid", "_density")
    
    # Define available options for probability density functions
    DISTRIBUTIONS = {'UNIFORM': UNIFORM, 'GAUSSIAN': GAUSSIAN, 'GMIXTURE': GMIXTURE,}
    
    def __init__(self, Ms, ab_max, ab_res, dH, sat_tol=1e-3, distribution="GAUSSIAN", dist_params=None, units="SI"):
        super().__init__(dH, units)

        # Raise exceptions for invalid input parameters
        if Ms<=0: raise ValueError("\"Ms\" (saturation magnetization) must be a positive, non-zero number")
        if ab_max<=0: raise ValueError("\"ab_max\" (maximum hysteron grid value) must be a positive, non-zero number")
        if ab_res<=0: raise ValueError("\"ab_res\" (hysteron grid resolution) must be a positive, non-zero number")
        if sat_tol<=0: raise ValueError("\"sat_tol\" (saturation tolerance) must be a positive, non-zero number")

        # Set valid input parameters
        self._Ms = Ms
        self._ab_max = ab_max
            
        # Construct the scaled Preisach grid as a flat column of alpha/beta pairs
        n_grid = int(round(ceil(2*self._ab_max/ab_res)))
        self._ab_grid = -1+2*array([[a,b] for a in range(n_grid+1) for b in range(n_grid+1)])/n_grid

        # Raise exceptions for invalid input distributions & parameters
        if isinstance(distribution, str):
            if distribution.upper() not in self.DISTRIBUTIONS:
                raise ValueError("Invalid choice of distribution, must be one of: {:s}".format(', '.join(self.DISTRIBUTIONS)))
            else:
                try:
                    self._density = self.DISTRIBUTIONS[distribution.upper()](self._ab_grid, *dist_params)
                except:
                    raise RuntimeError("Unable to call {:s} with arguments in \"dist_params\"".format(distribution.upper()))
        elif isinstance(distribution, ndarray):
            if distribution.min()<0:
                raise ValueError("Custom distributions must contain only non-negative values")
            elif distribution.shape[0]!=self._ab_grid.shape[0]:
                raise ValueError("Custom distributions must be the same shape as the Preisach grid")
            else:
                self._density = distribution
        else:
            raise ValueError("\"distribution\" must be a str or a numpy.ndarray")        

        # Reduce the Preisach density to its upper half, normalize, & rescale
        self._density[self._ab_grid[:,0]>=self._ab_grid[:,1]] = 0        
        self._density *= self._Ms/self._density.sum()
        self._ab_grid *= self._ab_max
        
        # Compute the major H-M curve & its critical points
        self._get_major(sat_tol)
                
    def _convert_params(self, cf):
        """Converts field units of Preisach parameters"""
        
        self._Ms *= cf
        self._ab_max *= cf
        self._ab_grid *= cf
        self._density *= cf
        return 0
        
    def _get_grid(self, ab0):
        R = zeros(len(self._ab_grid))
        R[-self._ab_grid[:,0] >= self._ab_grid[:,1]-ab0] = 1
        R[-self._ab_grid[:,0] < self._ab_grid[:,1]-ab0] = -1
        return R
    
    def _get_major(self, sat_tol):
        """Computes the major H-M curve & its critical points"""
        
        # Increase applied field until saturation is reached (initial curve)
        H_init = [0]
        M_init = [1e-6]
        R = self._get_grid(1e-6)
        while H_init[-1]<=self._ab_max:
            H_init.append(H_init[-1] + self.dH)
            R[H_init[-1]>self._ab_grid[:, 1]] = 1
            M_init.append(self._density@R)
            
        n_curve = int(ceil(2*max(H_init)/self.dH))+1
        
        # Decrease applied field until negative saturation is reached
        H_upper = []
        M_upper = []
        for i in range(n_curve):
            H_upper.append(H_init[-1]-self.dH*i)
            R[H_upper[-1]<self._ab_grid[:, 0]] = -1
            M_upper.append(self._density@R)
            
        # Increase applied field until saturation is reached again
        H_lower = []
        M_lower = []
        for i in range(n_curve):
            H_lower.append(H_upper[-1]+self.dH*i)   
            R[H_lower[-1]>self._ab_grid[:, 1]] = 1
            M_lower.append(self._density@R)
              
        # Combine H-M curves & determine critical points
        self.H_major = [h for h in H_init+H_upper+H_lower]
        self.M_major = [m for m in M_init+M_upper+M_lower]
        self.remanence = [r for r in [M_upper[H_upper.index(0.)], M_lower[H_lower.index(0.)]]]
        cuID = int(diff(sign(M_lower)).nonzero()[0])
        clID = int(diff(sign(M_upper)).nonzero()[0])
        self.coercivity = [c for c in [sum(H_lower[cuID:cuID+2])/2, sum(H_upper[clID:clID+2])/2]]
        
    def path(self, H_path, M0):
        """Computes hysteresis curves along applied field paths"""
        
        H = []
        M = []
        R = self._get_grid(self._ab_max*(M0/max(self.M_major)))
        
        # Loop over magnetizing field paths
        for Hp in H_path:
            deltaH = Hp[-1]-Hp[0]
            dH = self.dH*sign(deltaH)
            Np = int(ceil(abs(deltaH/dH)))+1
            
            # Compute the hysteresis along this path
            for n in range(Np):
                H.append(Hp[0] + n*dH)
                if dH>0: R[H[-1]>=self._ab_grid[:, 1]] = 1
                else: R[H[-1]<self._ab_grid[:, 0]] = -1
                M.append(self._density@R)
        
        return H, M

