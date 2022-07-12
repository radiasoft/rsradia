### hysteresis_models.py
### Morgan Henderson, July 2022
### A collection of magnetic hysteresis models

from numpy import pi, sinh, cosh, sin, cos, exp, array, zeros
from integrators import KRON15, EULER, RK4, RK45
from density_functions import gaussian

# Define a minimum magnetizing field strength resolution
delta_H = 1e-3

# Define the vacuum permeability constant
mu0 = 4*pi*1e-7

# Define numerical integrator options
INTEGRATORS = {
    'EULER': EULER,
    'RK4': RK4,
    'RK45': RK45,
}

# The jiles_atherton class requires an input parameter dictionary with keys:
#
# ja_params = {
#    alpha: Domain coupling strength
#    a: Domain wall density
#    Ms: Saturation magnetization of material
#    k: Pinning site breaking energy
#    c: Magnetization reversability
#    wa: Relative weight of anisotropic contributions to magnetization
#    Ka: Average anisotropy energy density
#    psi: Angle between anisotropy easy axis & magnetizing field
# }
#
# The default parameters are those for grade 20 steel as given by 
# Podbereznaya et al (https://doi.org/10.3103/S1068371219010115).
    
# A Jiles-Atherton hysteresis model
class jiles_atherton:
    
    # Initializes the model
    def __init__(self, **kwargs):
        
        # Define default parameters for the model
        ja_params = {
            'alpha': 4.93e-4,
            'a': 399,
            'Ms': 1.35e6,
            'k': 300,
            'c': 0.120,
            'wa': 0,
            'Ka': 0,
            'psi': 0,
        }
        
        # Assign input parameters to the model or use defaults
        for key in ja_params: setattr(self, key,kwargs.get(key, ja_params[key]))
        self.Ka /= mu0
                                      
        # Enforce isotropy if specified by EITHER Ka or wa
        if (not self.wa) or (not self.Ka): self.Ka = self.psi = self.wa = 0
        
        # Define anhysteretic magnetizations & derivatives as lambda functions
        self.Mai = lambda He: self.Ms*(cosh(He/self.a)/sinh(He/self.a)-self.a/He)
        self.dMai = lambda He: self.Ms*(self.a/He**2-1/(self.a*sinh(He/self.a)**2))
        self.Maa = lambda He: self.Ms*self.aniso_ints(He)
        self.dMaa = lambda He: (self.Maa(He+delta_H)-self.Maa(He-delta_H))/2e-3

    # Computes & takes the ratio of integrals for anisotropic, anhysteretic magnetization (Maa)
    def aniso_ints(self, He):
        
        # Define lambda functions for E terms
        E1 = lambda theta : (He*cos(theta)-(self.Ka/self.Ms)*sin(self.psi-theta)**2)/self.a
        E2 = lambda theta : (He*cos(theta)-(self.Ka/self.Ms)*sin(self.psi+theta)**2)/self.a
        
        # Define lambda functions for the integrands in the numerator & denominator
        f1 = lambda theta : exp((E1(theta)+E2(theta))/2)*sin(theta)
        f2 = lambda theta : f1(theta)*cos(theta)
        
        # Compute the integrals & return their ratio
        return KRON15(f2, (0, pi))/KRON15(f1, (0, pi))
             
    # The Jiles-Atherton differential equation for magnetization (returns dM/dH)
    def ja_diffeq(self, H, M, dH_sign):
        
        # Compute the effective magnetic field
        He = H + self.alpha*M
        
        # Compute & add anhysteretic, isotropic magnetization contributions
        Ma = (1-self.wa)*self.Mai(He)
        dMa = (1-self.wa)*self.dMai(He)
        
        # Compute & add anisotropic contributions (if any)
        if self.wa:
            Ma += self.wa*self.Maa(He)
            dMa += self.wa*self.dMaa(He)
            
        # Compute & return the magnetization differential
        return ((Ma-M)/(dH_sign*self.k-self.alpha*(Ma-M))+self.c*dMa)/(1+self.c)
        
    # Computes & returns magnetic hysteresis curves
    def get_hysteresis(self, **kwargs):
        
        # Assign hysteresis calculation options or use defaults
        hyst_opts = {
            'integrator': 'RK4',
            'H_lim': 1e3,
            'dH': 0,
        }
        for key in hyst_opts: hyst_opts[key] = kwargs.get(key,hyst_opts[key])
        if hyst_opts['dH']==0: hyst_opts['dH'] = 2*hyst_opts['H_lim']/1e3
        if hyst_opts['integrator'] not in INTEGRATORS.keys():
            print("WARNING: Integrator choice not understood, using RK4.")
            INT = RK4
        else: INT = INTEGRATORS[hyst_opts['integrator']]            
        H_lim, dH = hyst_opts['H_lim'], hyst_opts['dH']
            
        # Initialize hysteresis variables
        N = int(H_lim/dH)
        H = zeros(5*N)
        M = zeros(5*N)
        
        # Compute the hysteresis during initial magnetization
        H[:N] = array(range(N))*dH
        M[:N] = INT(0,(delta_H,H_lim), dH, self.ja_diffeq,1)
        
        # Compute the hysteresis during reverse magnetization
        H[N:3*N] = H_lim-array(range(2*N))*dH
        M[N:3*N] = INT(M[N-1], (H_lim-dH,-H_lim), -dH, self.ja_diffeq,-1)
        
        # Compute the hysteresis during re-magnetization
        H[3*N:] = -H_lim+array(range(2*N))*dH
        M[3*N:] = INT(M[3*N-1], (-H_lim+dH,H_lim), dH, self.ja_diffeq,1)
        
        return H, mu0*M
    
    def fit(self, H_true, B_true):
        None

# The preisach class requires an input parameter dictionary with keys:
#
# pr_params = {
#    Ms: Saturation magnetization of the magnet
#    grid: Grid of alpha/beta values 
#    density: Preisach density for alpha/beta values
# }

# A Preisach hysteresis model
class preisach:
    def __init__(self, **kwargs):
        
        # Define default parameters for the model
        pr_params = {
            'Ms': 1.3e6,
            'grid' : 1e3*(-1+2*array([[a,b] for a in range(100) for b in range(100)])/100),
            'density': array([1,-1]),
        }
        
        # Assign input parameters to the model or use defaults
        for key in pr_params: setattr(self,key,kwargs.get(key,pr_params[key]))
        if self.density.sum()==0: self.density = gaussian(self.grid,array([0,0]),array([[.1,0],[0,.1]]))
        self.density[self.grid[:,0]>=self.grid[:,1]] = 0        
        self.density /= self.density.sum()
        
    '''
    def hyst_grid(self): None
    def hyst_sample(self): None
    '''
            
    def get_hysteresis(self, **kwargs):
        
        # Assign hysteresis calculation options or use defaults
        hyst_opts = {
            'H_lim': 1e3,
            'dH': 0,
        }
        for key in hyst_opts: hyst_opts[key] = kwargs.get(key,hyst_opts[key])
        if hyst_opts['dH']==0: hyst_opts['dH'] = 2*hyst_opts['H_lim']/1e3           
        H_lim, dH = hyst_opts['H_lim'], hyst_opts['dH']
        
        # Initialize hysteresis variables
        N = int(H_lim/dH)
        H = zeros(4*N)
        M = zeros(4*N)+1
        R = zeros(self.density.shape)+1
        
        # Compute the hysteresis
        for h in range(2*N):
            H[h] = H_lim-(h+1)*dH
            R[H[h]<=self.grid[:,1]] = 0
            M[h] = (self.density*R).sum()
        for h in range(2*N):
            H[2*N+h] = -H_lim+(h+1)*dH
            R[H[2*N+h]>=self.grid[:,0]] = 1
            M[2*N+h] = (self.density*R).sum()
            
        return H, mu0*2*self.Ms*(M-.5)
    
    #
    def fit(self, H_true, B_true):
        None