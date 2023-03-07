### jiles_atherton.py
### February 2023

from numpy import sinh, cosh, sin, cos, exp, diff, sign

from .. import PI, MU0
from ..utils.integrators import *
from . import KRON_SPLITS

from .hysteresis_model import HysteresisModel

class JilesAtherton(HysteresisModel):
    """
    A Jiles-Atherton model for magnetic hysteresis

    parameters:
        alpha: Domain coupling strength
        a: Domain wall density (A/m or T)
        Ms: Saturation magnetization of material (A/m or T)
        k: Pinning site breaking energy (A/m or T)
        c: Magnetization reversability
        wa: Relative weight of anisotropic contributions to magnetization
        Ka: Average anisotropy energy density (J/m^3)
        psi: Offset angle between anisotropy easy axis & magnetizing field
    """
    
    # Define options for numerical integrators used by the model
    INTEGRATORS = {'EULER': EULER,'RK4': RK4,'RK45': RK45,}
    
    # List slotted members of JilesAtherton (non-extendable)
    __slots__ = ("_alpha", "_a", "_Ms", "_k", "_c", "_wa", "_Ka", "_psi", "_integrate")
    
    def __init__(self, alpha, a, Ms, k, c, dH, wa=0, Ka=0, psi=0, sat_tol=1e-3, integrator="RK4", units="SI"):
        super().__init__(dH, units)

        # Raise exceptions for invalid input parameters
        if alpha<=0: raise ValueError("\"alpha\" (domain coupling strength) must be a positive, non-zero number")
        if a<=0: raise ValueError("\"a\" (domain wall density) must be a positive, non-zero number")
        if Ms<=0: raise ValueError("\"Ms\" (saturation magnetization) must be a positive, non-zero number")
        if k<=0: raise ValueError("\"k\" (pin-breaking energy) must be a positive, non-zero number")
        if c<=0: raise ValueError("\"c\" (magnetization reversability) must be a positive, non-zero number")
        if wa<0 or wa>1: raise ValueError("\"wa\" (anisotropy weight) must be at least zero and at most one")
        if Ka<0: raise ValueError("\"Ka\" (average anisotropy density) must be zero or a positive number")
        if psi<0 or psi>PI: raise ValueError("\"psi\" (anisotropy easy axis offset angle) must be at least zero and at most pi")
        if integrator.upper() not in self.INTEGRATORS:
            raise ValueError("\"integrator\" must be one of {:s}".format(', '.join(self.INTEGRATORS)))

        # Set input parameters
        self._alpha = alpha
        self._a = a/self.units
        self._Ms = Ms/self.units
        self._k = k/self.units
        self._c = c
        if wa and Ka:
            self._wa = wa
            self._Ka = (Ka/MU0)
            self._psi = psi
        else: 
            self._Ka = self._psi = self._wa = 0
        self._integrate = self.INTEGRATORS[integrator.upper()]
        
        # Compute the major H-M curve & its critical points
        self._get_major(sat_tol)
    
    def _get_aniso(self, He):
        """Computes the anhysteretic, anisotropic magnetization"""
        
        def _den_int(self, theta):
            """Computes the integrand of the denominator in the anhysteretic, anisotropic magnetization"""
            
            E1 = (He*cos(theta)-(self._Ka/(self._Ms*self._a))*sin(self._psi-theta)**2)
            E2 = (He*cos(theta)-(self._Ka/(self._Ms*self._a))*sin(self._psi+theta)**2)
            return exp((E1+E2)/2)*sin(theta)
    
        def _num_int(self, theta):
            """Computes the integrand of the numerator in the anhysteretic, anisotropic magnetization"""

            return self._den_int(theta)*cos(theta)
        
        # Compute the 
        Maa = self._Ms*KRON15(_num_int, (0, PI), KRON_SPLITS)/KRON15(_den_int, (0, PI), KRON_SPLITS)
        dMaa = (Maa(He+self.dH)-Maa(He-self.dH))/(2*self.dH)
        
        return Maa, dMaa
             
    def _diffeq(self, H, M, delta):
        """Computes the Jiles-Atherton differential equation for magnetization (returns dM/dH)"""
        
        # Compute the effective magnetic field
        He = H + self._alpha*M
        if He==0: He+=1e-6*delta
        
        # Compute & add anhysteretic, isotropic magnetization contributions
        Ma = (1-self._wa)*self._Ms*(cosh(He/self._a)/sinh(He/self._a)-self._a/He)
        dMa = (1-self._wa)*self._Ms*(self._a/He**2-1/(self._a*sinh(He/self._a)**2))
        
        # Compute & add anfisotropic contributions (if any)
        if self._wa:
            Maa, dMaa = self._get_aniso(He)
            Ma += self._wa*Maa
            dMa += self._wa*dMaa
            
        # Compute & return the magnetization differential
        return ((Ma-M)/(delta*self._k-self._alpha*(Ma-M))+self._c*dMa)/(1+self._c)
    
    def _get_major(self, sat_tol):
        """Computes the major H-M curve & its critical points"""
        
        # Increase applied field until saturation is reached (initial curve)
        delta_M = 1
        H_init = [0]
        M_init = [1e-6]
        while (delta_M > sat_tol):
            H_init.append(H_init[-1]+self.dH)            
            M_init.append(self._integrate(M_init[-1], H_init[-2:], self.dH, self._diffeq, 1)[-1])
            delta_M = abs((M_init[-1]-M_init[-2])/M_init[-2])
                   
        # Cycle to negative saturation & return to positive saturation (upper & lower curves)
        H_upper = [H_init[-1]-self.dH*n for n in range(2*len(H_init)+1)]
        M_upper = self._integrate(M_init[-1], [H_upper[0], H_upper[-1]], -self.dH, self._diffeq, -1)      
        H_lower = [H_upper[-1]+self.dH*n for n in range(2*len(H_init)+1)]
        M_lower = self._integrate(M_upper[-1], [H_lower[0], H_lower[-1]], self.dH, self._diffeq, 1)   
        
        # Combine H-M curves & determine critical points
        self.H_major = [h*self.units for h in H_init+H_upper+H_lower]
        self.M_major = [m*self.units for m in M_init+M_upper+M_lower]
        self.remanence = [r*self.units for r in [M_upper[H_upper.index(0.)], M_lower[H_lower.index(0.)]]]
        cuID = int(diff(sign(M_lower)).nonzero()[0])
        clID = int(diff(sign(M_upper)).nonzero()[0])
        self.coercivity = [c*self.units for c in [sum(H_lower[cuID:cuID+2])/2, sum(H_upper[clID:clID+2])/2]]
        
    def path(self, H_path, M0):
        """Computes a hysteresis curve along a path of applied fields"""
        
        H = []
        M = []
        M0 /= self.units
        
        # Loop over magnetizing field paths
        for p in range(len(H_path)):
            
            # Define starting & stopping field strengths, field step, & number of steps
            Hp = [h/self.units for h in H_path[p]]
            dH = self.dH/self.units if Hp[-1] > Hp[0] else -self.dH/self.units
            n_path = int(round(ceil((Hp[-1]-Hp[0])/dH)))
            
            # Compute hysteresis values along the path
            H += [Hp[0]+dH*n for n in range(n_path+1)]
            M += self._integrate(M0, [Hp[0], H[-1]], dH, self._diffeq, sign(dH)) 
            M0 = M[-1]
        
        return [h*self.units for h in H], [m*self.units for m in M]
    