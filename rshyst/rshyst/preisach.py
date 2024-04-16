"""preisach.py

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

from rshyst import distributions, HysteresisModel

class Preisach(HysteresisModel):
    """A Preisach model of magnetic hysteresis"""
    
    def __init__(
            self, Ms, ab_max, ab_res, dH, 
            sat_tol=1e-3, distribution="GAUSSIAN", dist_params=None, units="SI"
    ):
        """Args:
          * Ms: Saturation magnetization of the magnet
          * ab_max: Maximum magnitude of hysteron on/off field values
          * ab_res: Resolution of alpha/beta grid (point spacing)
          * distribution: Type of weight distribution used
          * dist_params: Parameters of the weight distribution
        """

        super().__init__(dH, units)

        # Set valid input parameters
        self._Ms = Ms
        self._ab_max = ab_max
            
        # Construct the scaled Preisach grid as a flat column of alpha/beta pairs
        n_ab = int(round(np.ceil(2 * self._ab_max / ab_res)))
        self._ab = -1 + 2 * np.array([[a, b] for a in range(n_ab + 1) for b in range(n_ab + 1)]) / n_ab

        # Raise exceptions for invalid input distributions & parameters
        distribution = getattr(distributions, distribution)
        self._density = distribution(self._ab, *dist_params)     

        # Reduce the Preisach density to its upper half, normalize, & rescale
        self._density[self._ab[:,0]>=self._ab[:,1]] = 0        
        self._density *= self._Ms / self._density.sum()
        self._ab *= self._ab_max
        
        # Compute the major H-M curve & its critical points
        self._get_major(sat_tol)
                
    def _convert_params(self, cf):
        """Converts field units of Preisach parameters"""
        
        self._Ms *= cf
        self._ab_max *= cf
        self._ab *= cf
        self._density *= cf
        
        return 0
        
    def _get_grid(self, ab0):

        R = np.zeros(len(self._ab))
        R[-self._ab[:,0] >= self._ab[:,1]-ab0] = 1
        R[-self._ab[:,0] < self._ab[:,1]-ab0] = -1
        return R
    
    def _get_major(self, sat_tol):
        """Computes the major H-M curve & its critical points"""
        
        # Increase applied field until saturation is reached (initial curve)
        H_init = [0]
        M_init = [1e-6]
        R = self._get_grid(1e-6)
        while H_init[-1]<=self._ab_max:
            H_init.append(H_init[-1] + self.dH)
            R[H_init[-1]>self._ab[:, 1]] = 1
            M_init.append(self._density @ R)
            
        n_curve = int(np.ceil(2*max(H_init)/self.dH))+1
        
        # Decrease applied field until negative saturation is reached
        H_upper = []
        M_upper = []
        for i in range(n_curve):
            H_upper.append(H_init[-1]-self.dH*i)
            R[H_upper[-1]<self._ab[:, 0]] = -1
            M_upper.append(self._density@R)
            
        # Increase applied field until saturation is reached again
        H_lower = []
        M_lower = []
        for i in range(n_curve):
            H_lower.append(H_upper[-1]+self.dH*i)   
            R[H_lower[-1]>self._ab[:, 1]] = 1
            M_lower.append(self._density@R)
              
        # Combine H-M curves & determine critical points
        self.H_major = [h for h in H_init+H_upper+H_lower]
        self.M_major = [m for m in M_init+M_upper+M_lower]
        self.remanence = [r for r in [M_upper[H_upper.index(0.)], M_lower[H_lower.index(0.)]]]
        cuID = int(np.diff(np.sign(M_lower)).nonzero()[0])
        clID = int(np.diff(np.sign(M_upper)).nonzero()[0])
        self.coercivity = [c for c in [sum(H_lower[cuID:cuID+2])/2, sum(H_upper[clID:clID+2])/2]]
        
    def path(self, H_path, M0):
        """Computes hysteresis curves along applied field paths"""
        
        H = []
        M = []
        R = self._get_grid(self._ab_max * (M0 / self._Ms))
        
        # Loop over magnetizing field paths
        for Hp in H_path:
            deltaH = Hp[-1]-Hp[0]
            dH = self.dH*np.sign(deltaH)
            Np = int(np.ceil(abs(deltaH/dH)))+1
            
            # Compute the hysteresis along this path
            for n in range(Np):
                H.append(Hp[0] + n*dH)
                if dH>0: R[H[-1]>=self._ab[:, 1]] = 1
                else: R[H[-1]<self._ab[:, 0]] = -1
                M.append(self._density@R)
        
        return H, M
