"""model.py

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

import pickle
import numpy as np

class HysteresisModel:
    """Base class for magnetic hysteresis models"""    

    def __init__(self, dH):
            
        # Set input parameters
        self.dH = dH
        
    def _get_major(self, sat_tol):
        """VIRTUAL: computes the major H-M curve & its critical points"""

        raise NotImplemented("This is a virtual method. Stop.")
        
    def _convert_params(self, cf):
        """Converts model-specific parameter field units"""

        return 0
        
    def path(self, H_path, M0):
        """VIRTUAL: computes hysteresis curves along applied field paths"""

        raise NotImplemented("This is a virtual method. Stop.")
                
    def point(self, H_point, curve="upper"):
        """Computes magnetization on a curve of the major loop at a particular field strength"""
        
        # Set indices based on chosen curve
        if curve.lower()=='initial':
            indices = (0, int(.2*len(self.H_major)))
        elif curve.lower()=='upper':
            indices = (int(.2*len(self.H_major)), int(.6*len(self.H_major)))
        elif curve.lower()=='lower':
            indices = (int(.6*len(self.H_major)), int(len(self.H_major)))
        else:
            raise ValueError("\"curve\" must be one of: initial, upper, lower")
        
        # Select points along curve
        H = self.H_major[indices[0]:indices[1]]
        M = self.M_major[indices[0]:indices[1]]
        
        # Find curve point closest to prescribed point
        H_dists = [abs(h - H_point) for h in H]
        neighbors = np.argsort(H_dists)[:2]
            
        # Interpolate magnetic flux density between neighboring points
        return sum([M[pt]*H_dists[pt] for pt in neighbors])/self.dH

    def save(self, path="./ja_model.pkl"):
        """Saves a hysteresis model object as a .pkl file"""

        with open(path,'wb') as file:
            pickle.dump(self, file)
    
    @classmethod
    def load(cls, path):
        """Loads a hysteresis model object from a pkl file"""
        
        with open(path, 'rb') as file:
            return pickle.load(file)
        