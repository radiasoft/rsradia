### jiles_atherton.py
### November 2022
### A Jiles-Atherton model for the hysteresis of a magnetic material

from numpy import array, zeros, sinh, cosh, sin, cos, exp, ceil, argsort
from pickle import dump, load

from .. import PI, MU0
from .integrators import *

# Define numerical integrator options
INTEGRATORS = {
    'EULER': EULER,
    'RK4': RK4,
    'RK45': RK45,
}

class JilesAtherton:
    '''
    A Jiles-Atherton model of magnetic hysteresis

    parameters:
        alpha: Domain coupling strength
        a: Domain wall density
        Ms: Saturation magnetization of material
        k: Pinning site breaking energy
        c: Magnetization reversability
        wa: Relative weight of anisotropic contributions to magnetization
        Ka: Average anisotropy energy density
        psi: Offset angle between anisotropy easy axis & magnetizing field
    '''
    
    def __init__(self, alpha, a, Ms, k, c, wa=0, Ka=0, psi=0, dH=1, integrator='RK4'):

        # Raise exceptions for invalid input parameters
        if alpha<=0: raise ValueError("\"alpha\" (domain coupling strength) must be a positive, non-zero number")
        if a<=0: raise ValueError("\"a\" (domain wall density) must be a positive, non-zero number")
        if Ms<=0: raise ValueError("\"Ms\" (saturation magnetization) must be a positive, non-zero number")
        if k<=0: raise ValueError("\"k\" (pin-breaking energy) must be a positive, non-zero number")
        if c<=0: raise ValueError("\"c\" (magnetization reversability) must be a positive, non-zero number")
        if wa<0 or wa>1: raise ValueError("\"wa\" (anisotropy weight) must be at least zero and at most one")
        if Ka<0: raise ValueError("\"Ka\" (average anisotropy density) must be zero or a positive number")
        if psi<0 or psi>PI: raise ValueError("\"psi\" (anisotropy easy axis offset angle) must be at least zero and at most pi")
        if dH<=0: raise ValueError("\"dH\" (field step-size) must be a positive, non-zero number")
        if integrator.upper() not in INTEGRATORS.keys():
            raise ValueError("\"integrator\" must be one of {:s}".format(', '.join(INTEGRATORS.keys())))

        # Set valid input parameters
        self._alpha = alpha
        self._a = a
        self._Ms = Ms
        self._k = k
        self._c = c
        self._wa = wa
        self._Ka = Ka
        self._psi = psi
        self._dH = dH
        self._integrator= INTEGRATORS[integrator.upper()]      

        # Enforce isotropy if specified by EITHER Ka or wa parameters
        self._Ka = Ka/MU0
        if (not self._wa) or (not self._Ka): self._Ka = self._psi = self._wa = 0        
        
        # Compute the major hysteresis curve for the model (sets self.H_major, self.B_major)
        self._get_major()
        
        # Define the remanant magnetic field values for the model
        rem_points = [int(.4*len(self.H_major)), int(.8*len(self.H_major))]
        self.B_rem = array([self.B_major[rpt] for rpt in rem_points])
            
    # Computes the anhysteretic, isotropic magnetization
    def _Mai(self, He):
        return self._Ms*(cosh(He/self._a)/sinh(He/self._a)-self._a/He)
    
    # Computes the anhysteretic, isotropic magnetization derivative
    def _dMai(self, He):
        return self._Ms*(self._a/He**2-1/(self._a*sinh(He/self._a)**2))
    
    # Computes the anhysteretic, anisotropic magnetization
    def _Maa(self, He):
        
        # Define the integrand of the denominator
        def _denominator_integrand(theta):
            E1 = (He*cos(theta)-(self._Ka/self._Ms)*sin(self._psi-theta)**2)/self._a
            E2 = (He*cos(theta)-(self._Ka/self._Ms)*sin(self._psi+theta)**2)/self._a
            return exp((E1+E2)/2)*sin(theta)
        
        # Define the integrand of the numerator
        def _numerator_integrand(theta):
            return _denominator_integrand(theta)*cos(theta)
        
        # Compute the integrals & return the result
        return self._Ms*KRON15(_numerator_integrand, (0, PI))/KRON15(_denominator_integrand, (0, PI))
    
    # Computes the anhysteretic, anisotropic mangetization derivative
    def _dMaa(self, He):
        return (self._Maa(He+self._dH)-self._Maa(He-self._dH))/(2*self._dH)
             
    # Computes the Jiles-Atherton differential equation for magnetization (returns dM/dH)
    def _diffeq(self, H, M, dH_sign):
        
        # Compute the effective magnetic field
        He = H + self._alpha*M
        
        # Compute & add anhysteretic, isotropic magnetization contributions
        Ma = (1-self._wa)*self._Mai(He)
        dMa = (1-self._wa)*self._dMai(He)
        
        # Compute & add anfisotropic contributions (if any)
        if self._wa:
            Ma += self._wa*self._Maa(He)
            dMa += self._wa*self._dMaa(He)
            
        # Compute & return the magnetization differential
        return ((Ma-M)/(dH_sign*self._k-self._alpha*(Ma-M))+self._c*dMa)/(1+self._c)
        
    # Computes the major magnetic hysteresis curve of the model
    def _get_major(self):
        
        # Initialize the applied field & magnetization
        H = [0]
        M = [1e-6]
        
        # Increase applied field until saturation is reached
        delta_M = 1
        while (delta_M > 2.5e-5):
            
            # Increase applied field & compute the magnetization 
            H.append(H[-1] + self._dH)
            M.append(self._integrator(M[-1], (H[-2], H[-1]), self._dH, self._diffeq, 1)[1])

            # Determine the relative change in magnetization over this step
            dMdH = abs(M[-1]-M[-2])/self._dH
            delta_M = dMdH/M[-2]
            
        n_curve = int(round(ceil(2*max(H)/self._dH)))
            
        # Decrease applied field until negative saturation is reached
        H += [max(H)-self._dH*n for n in range(n_curve+1)]
        M += list(self._integrator(M[-1], (max(H), -max(H)), -self._dH, self._diffeq, -1))
        
        # Increase applied field until saturation is reached again
        H += [min(H)+self._dH*n for n in range(n_curve+1)]
        M += list(self._integrator(M[-1], (min(H), max(H)), self._dH, self._diffeq, 1))        
            
        # Assign the hysteresis variables to model class members
        self.H_major = array(H)
        self.B_major = MU0*array(M)
    
    # Use the model to compute hysteresis curves along a path of applied fields
    def path(self, H_path, M0):
        
        # Determine numbers of hysteresis points & initialize hysteresis variables
        N_path = (ceil(H_path.ptp(axis=1))/self._dH).round().astype('int32')+1
        H = zeros(N_path.sum())
        M = zeros(N_path.sum())
        
        # Loop over magnetizing field paths
        for p in range(len(H_path)):
            
            # Define starting & stopping field strengths for this path
            Hstart, Hstop = H_path[p]
            sgn = 1 if Hstop > Hstart else -1
            
            # Determine starting index & number of points on this path
            N0 = N_path[:p].sum()
            N = N_path[p]
            
            # Compute the hysteresis along this path
            H[N0:N0+N] = Hstart + array(range(N))*self._dH*sgn
            M[N0:N0+N] = self._integrator(M0, (Hstart, Hstop), self._dH*sgn, self._diffeq, sgn)
            M0 = M[N0+N-1]
        
        return H, MU0*M
    
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
    
    # Saves a Jiles-Atherton model to a file
    def save(self, path="./ja_model.pkl"):
        with open(path,'wb') as file:
            dump(self, file)
    
    # Loads a Jiles-Atherton model from a file (possibly on instantiation)
    @classmethod
    def load(cls, path="./ja_model.pkl"):
        with open(path, 'rb') as file:
            return load(file)