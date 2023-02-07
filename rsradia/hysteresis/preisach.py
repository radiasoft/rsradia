### preisach.py
### November 2022
### A Preisach model for the hysteresis of a magnetic material

from numpy import array, zeros, ceil, argsort, ndarray, diff, sign
from pickle import dump, load

from . import MU0
from .distributions import *

# Define available options for probability density functions
DISTRIBUTIONS = {
    'UNIFORM': UNIFORM,
    'GAUSSIAN': GAUSSIAN,
    'GMIXTURE': GMIXTURE,
}

class Preisach:
    '''
    A Preisach model of magnetic hysteresis

    Parameters:
        Ms: Saturation magnetization of the magnet
        ab_max: Maximum magnitude of hysteron on/off field values
        ab_res: Resolution of alpha/beta grid (point spacing)
        distribution: Type of weight distribution used
        dist_params: Parameters of the weight distribution
    '''
    
    # Initializes a Preisach hysteresis model
    def __init__(self, Ms, ab_max, ab_res, dH=1, distribution="GAUSSIAN", dist_params=None):

        # Raise exceptions for invalid input parameters
        if Ms<=0: raise ValueError("\"Ms\" (saturation magnetization) must be a positive, non-zero number")
        if ab_max<=0: raise ValueError("\"ab_max\" (maximum hysteron grid value) must be a positive, non-zero number")
        if ab_res<=0: raise ValueError("\"ab_res\" (hysteron grid resolution) must be a positive, non-zero number")
        if dH<=0: raise ValueError("\"dH\" (field step-size) must be a positive, non-zero number")

        # Set valid input parameters
        self._Ms = Ms
        self._ab_max = ab_max
        self._ab_res = ab_res
        self._dH = dH
            
        # Construct the scaled Preisach grid
        n_grid = int(round(ceil(2*self._ab_max/self._ab_res)))
        self._grid = -1+2*array([[a,b] for a in range(n_grid+1) for b in range(n_grid+1)])/n_grid

        # Raise exceptions for invalid input distributions & parameters
        if isinstance(distribution, str):
            if distribution.upper() not in DISTRIBUTIONS.keys():
                raise ValueError("Invalid choice of distribution, must be one of: {:s}".format(', '.join(DISTRIBUTIONS.keys())))
            else:
                try:
                    self._density = DISTRIBUTIONS[distribution.upper()](self._grid, *dist_params)
                except:
                    raise RuntimeError("Unable to call {:s} with arguments in \"dist_params\"".format(distribution.upper()))
        elif isinstance(distribution, ndarray):
            if distribution.min()<0:
                raise ValueError("Custom distributions must contain only non-negative values")
            elif distribution.shape[0]!=self._grid.shape[0]:
                raise ValueError("Custom distributions must be the same shape as the Preisach grid")
            else:
                self._density = distribution
        else:
            raise ValueError("\"distribution\" must be a str or a numpy.ndarray")        

        # Reduce the Preisach density to its upper half & normalize, rescale the Preisach grid
        self._density[self._grid[:,0]>=self._grid[:,1]] = 0        
        self._density /= self._density.sum()
        self._grid *= self._ab_max
        
        # Compute the major hysteresis curve for the model (sets self.H_major, self.B_major)
        self._get_major()
        self._get_criticals()
            
    # Computes the initial & major 
    def _get_major(self):
        
        # Initialize the applied field & magnetization
        H = [0]
        M = [0]
        
        # Initialize hysteron values
        R = zeros(len(self._grid))
        R[-self._grid[:,0] >= self._grid[:,1]] = 1
        R[-self._grid[:,0] < self._grid[:,1]] = -1
        R[self._grid[:,0]>=self._grid[:,1]] = 0
        
        # Increase applied field until saturation is reached
        while (H[-1]<self._ab_max):
            H.append(H[-1] + self._dH)
            R[H[-1]>self._grid[:, 1]] = 1
            M.append(self._density@R)
            
        n_curve = int(round(ceil(2*max(H)/self._dH)))
        
        # Decrease applied field until negative saturation is reached
        for i in range(n_curve):
            H.append(H[-1] - self._dH)
            R[H[-1]<self._grid[:, 0]] = -1
            M.append(self._density@R)
            
        # Increase applied field until saturation is reached again
        for i in range(n_curve):
            H.append(H[-1] + self._dH)
            R[H[-1]>self._grid[:, 1]] = 1
            M.append(self._density@R)
              
        # Assign the hysteresis variables to model class members
        self.H_major = array(H)
        self.B_major = self._Ms*MU0*array(M)
        
    # Computes the material coercivity & remanent field
    def _get_criticals(self):        
        
        # Separate the major loop into upper & lower sections
        upper_inds = range(int(.2*len(self.H_major)), int(.6*len(self.H_major)))
        lower_inds = range(int(.6*len(self.H_major)), int(len(self.H_major)))
        H_upper = self.H_major[upper_inds]
        H_lower = self.H_major[lower_inds]
        B_upper = self.B_major[upper_inds]
        B_lower = self.B_major[lower_inds]
        
        # Determine the remanent field
        if 0 in H_upper:
            Bru = float(B_upper[H_upper==0])
        else:
            H0 = int(diff(sign(H_upper)).nonzero()[0])
            Bru = float(B_upper[H0:H0+2].sum()/2)
        if 0 in H_lower:
            Brl = float(B_lower[H_lower==0])
        else:
            H0 = int(diff(sign(H_lower)).nonzero()[0])
            Brl = float(B_lower[H0:H0+2].sum()/2)
        self.remanence = array([Brl, Bru])
            
        # Determine the coercivity
        if 0 in B_upper:
            Hcu = float(H_upper[B_upper==0])
        else:
            B0 = int(diff(sign(B_upper)).nonzero()[0])
            Hcu = float(H_upper[B0:B0+2].sum()/2)
        if 0 in B_lower:
            Hcl = float(H_lower[B_lower==0])
        else:
            B0 = int(diff(sign(B_lower)).nonzero()[0])
            Hcl = float(H_lower[B0:B0+2].sum()/2)
        self.coercivity = array([Hcu, Hcl])
        
    # Use the model to compute hysteresis curves along a path of applied fields
    def path(self, H_path, M0):
        
        # Determine numbers of hysteresis points & field change signs along path
        N_path = (ceil(H_path.ptp(axis=1))/self._dH).round().astype('int32')+1
        sgn_path = array([1 if Hs[1]>Hs[0] else -1 for Hs in H_path])
        
        # Initialize hysteresis variable containers
        H = zeros(N_path.sum())
        M = zeros(N_path.sum())
        
        # Initializehysteron values
        ab0 = (M0/self.B_major.max()/MU0)*self._ab_max
        R = zeros(len(self._grid))
        R[-self._grid[:,0] >= self._grid[:,1]-ab0] = 1
        R[-self._grid[:,0] < self._grid[:,1]-ab0] = -1
        
        # Loop over magnetizing field paths
        for p in range(len(H_path)):
            
            # Define starting & stopping field strengths for this path
            Hstart, Hstop = H_path[p]
            sgn = sgn_path[p]
            
            # Determine starting index & number of points on this path
            N0 = N_path[:p].sum()
            N = N_path[p]
            
            # Compute the hysteresis along this path
            for n in range(N):
                H[N0+n] = Hstart + n*self._dH*sgn
                if sgn>0: R[H[N0+n] > self._grid[:, 1]] = 1
                else: R[H[N0+n] < self._grid[:, 0]] = -1
                M[N0+n] = self._density@R
        
        return H, self._Ms*MU0*M
                
    # Use the model to compute magnetic flux density at a particular field strength
    def point(self, H_point, curve='upper'):
        
        # Indices for points along each curve
        curve_indices = {
            'initial': (0, int(.2*len(self.H_major))),
            'upper': (int(.2*len(self.H_major)), int(.6*len(self.H_major))),
            'lower': (int(.6*len(self.H_major)), int(len(self.H_major)))
        }
        
        # Select points along the desired curve
        if curve.lower() in curve_indices.keys():
            indices = curve_indices[curve.lower()]
        else:
            print("WARNING: Curve choice not understood, using upper curve.")
            indices = curve_indices['upper']
        H = self.H_major[indices[0]:indices[1]]
        B = self.B_major[indices[0]:indices[1]]
        
        # Find curve point closest to prescribed point
        H_dists = abs(H - H_point)
        neighbors = argsort(H_dists)[:2]
            
        # Interpolate magnetic flux density between neighboring points
        return sum([B[pt]*H_dists[pt] for pt in neighbors])/self._dH
    
    #
    @classmethod
    def fit(cls, **fit_params):
        """
        
        Parameters:
            
        """
        return 0
    
    # Saves a Preisach model to a file
    def save(self, path="./self.pkl"):
        with open(path,'wb') as file:
            dump(self, file)
    
    # Loads a Preisach model from a file (possibly on instantiation)
    @classmethod
    def load(cls, path="./self.pkl"):
        with open(path, 'rb') as file:
            return load(file)