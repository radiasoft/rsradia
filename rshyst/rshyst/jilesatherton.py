"""jilesatherton.py

LICENSE STATEMENT

Copyright 2024 RadiaSoft LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np

from rshyst.constants import MU0
from rshyst import integrators, HysteresisModel

class JilesAtherton(HysteresisModel):
    """A Jiles-Atherton magnetic hysteresis model"""
        
    def __init__(
        self, alpha, a, Ms, k, c, dH, wa=0, Ka=0, theta=0, phi=0,
        sat_tol=1e-3, integrator="RK4"
    ):
        """Args:
          * alpha: Domain coupling strength
          * a: Domain wall density (A/m or T)
          * Ms: Saturation magnetization of material (A/m or T)
          * k: Pinning site breaking energy (A/m or T)
          * c: Magnetization reversability
          * dH: Magnetizing field step-size ()
          * wa: Relative weight of anisotropic effects (default 0.)
          * Ka: Average anisotropy energy density (J/m^3, default 0)
          * theta: Easy axis polar angle (radians, default 0)
          * phi: Easy axis azimuthal angle (radians, default 0)
          * sat_tol: Largest relative change in magnetization at saturation (default 1e-3)
          * integrator: Integrator used for solving differential equation (default 'RK4')
        """

        super().__init__(dH)

        # Set isotropic input parameters
        self.alpha = alpha
        self.a = a
        self.Ms = Ms
        self.k = k
        self.c = c

        # Set anisotropic input parameters
        if wa and Ka:
            self.wa = wa
            self.Ka = (Ka/MU0)
            self.theta = theta
            self.phi = phi
        else: 
            self.Ka = self.wa = self.theta = self.phi = self.psi = 0

        # Determine angles between easy axis & cartesian axes
        easy_unit = np.array([
            np.sin(self.theta) * np.cos(self.phi),
            np.sin(self.theta) * np.sin(self.phi),
            np.cos(self.theta)
        ])
        self.psi = np.array([np.arccos(easy_unit @ u) for u in np.eye(3)])

        #
        self.integrator = getattr(integrators, integrator)
        
        # Compute the major H-M curve & its critical points
        self.get_major(sat_tol)

    def get_Maa(self, He):
        """Computes anisotropic magnetizaiton as a function of applied field"""
        
        def den_int(self, theta):
            """Computes the integrand of the denominator in the anhysteretic, anisotropic magnetization"""
            
            A = He * np.cos(theta) - self.Ka / (self.Ms * self.a)
            E1 = A * np.sin(self.psi - theta)**2
            E2 = A * np.sin(self.psi + theta)**2
            Id = np.exp((E1 + E2) / 2) * np.sin(theta)

            return Id
    
        def num_int(self, theta):
            """Computes the integrand of the numerator in the anhysteretic, anisotropic magnetization"""

            In = self.den_int(theta, self.psi) * np.cos(theta)

            return In
        
        # Compute the 
        phi_lims = (0, np.pi)
        num = integrators.KRON15(num_int, phi_lims, 10)
        den = integrators.KRON15(den_int, phi_lims, 10)
        Maa = self.Ms * num / den

        return Maa

    def _get_aniso(self, He):
        """Computes the anhysteretic, anisotropic magnetization and its local derivative"""

        Maa = self.get_Maa(He)
        dMaa = (self.get_Maa(He + self.dH) - self.get_Maa(He - self.dH)) / (2 * self.dH)
        
        return Maa, dMaa
             
    def diffeq(self, H, M, delta):
        """Computes the Jiles-Atherton differential equation for magnetization

        Args:
          * H: Current applied magnetic field
          * M: Current magnetization
          * delta: Sign of current change in applied magnetic field

        Returns:
          * dMdH: Magnetization differential
        """
        
        # Compute the effective magnetic field
        He = H + self.alpha*M
        if He==0: He+=1e-9*delta
        
        # Compute & add anhysteretic, isotropic magnetization contributions
        Ma = (1-self.wa)*self.Ms*(np.cosh(He/self.a)/np.sinh(He/self.a)-self.a/He)
        dMa = (1-self.wa)*self.Ms*(self.a/He**2-1/(self.a*np.sinh(He/self.a)**2))
        
        # Compute & add anfisotropic contributions (if any)
        if self.wa:
            Maa, dMaa = self.get_aniso(He)
            Ma += self.wa*Maa
            dMa += self.wa*dMaa
            
        # Compute & return the magnetization differential
        dMdH = ((Ma-M)/(delta*self.k-self.alpha*(Ma-M))+self.c*dMa)/(1+self.c)

        return dMdH
    
    def get_major(self, sat_tol):
        """Computes the major H-M curve & its critical points"""
        
        # Increase applied field until saturation is reached (initial curve)
        delta_M = 1
        H_init = [0]
        M_init = [1e-6]
        while (delta_M > sat_tol):
            H_init.append(H_init[-1]+self.dH)            
            M_init.append(self.integrator(M_init[-1], H_init[-2:], self.dH, self.diffeq, 1)[-1])
            delta_M = abs((M_init[-1]-M_init[-2])/M_init[-2])
                   
        # Cycle to negative saturation & return to positive saturation (upper & lower curves)
        H_upper = [H_init[-1]-self.dH*n for n in range(2*len(H_init)+1)]
        M_upper = self.integrator(M_init[-1], [H_upper[0], H_upper[-1]], -self.dH, self.diffeq, -1)      
        H_lower = [H_upper[-1]+self.dH*n for n in range(2*len(H_init)+1)]
        M_lower = self.integrator(M_upper[-1], [H_lower[0], H_lower[-1]], self.dH, self.diffeq, 1)   
        
        # Combine H-M curves & determine critical points
        self.H_major = np.array([h for h in H_init+H_upper+H_lower], dtype=float)
        self.M_major = np.array([m for m in M_init+M_upper+M_lower], dtype=float)
        self.remanence = np.array([r for r in [M_upper[H_upper.index(0.)], M_lower[H_lower.index(0.)]]], dtype=float)
        cuID = int(np.diff(np.sign(M_lower)).nonzero()[0])
        clID = int(np.diff(np.sign(M_upper)).nonzero()[0])
        self.coercivity = [c for c in [sum(H_lower[cuID:cuID+2])/2, sum(H_upper[clID:clID+2])/2]]
        
    def path(self, H_path, M0):
        """Computes a hysteresis curve along a path of applied fields"""
        
        H = []
        M = []
        
        # Loop over magnetizing field paths
        for p in range(len(H_path)):
            
            # Define starting & stopping field strengths, field step, & number of steps
            Hp = [h for h in H_path[p]]
            dH = self.dH if Hp[-1] > Hp[0] else -self.dH
            n_path = int(round(np.ceil((Hp[-1]-Hp[0])/dH)))
            
            # Compute hysteresis values along the path
            H += [Hp[0]+dH*n for n in range(n_path+1)]
            M += self.integrator(M0, [Hp[0], H[-1]], dH, self.diffeq, np.sign(dH)) 
            M0 = M[-1]
        
        return np.array(H, dtype=float), np.array(M, dtype=float)
