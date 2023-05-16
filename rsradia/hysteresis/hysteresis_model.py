### hysteresis_model.py
### February 2023

from .. import MU0

from numpy import argsort, diff, sign, ceil
from pickle import dump, load

class HysteresisModel:
    """A base class for models of magnetic hysteresis"""
    
    UNITS = {"SI":1, "MM":1e-3, "TESLA":MU0, "GAUSS":MU0*1e4}
    
    # Slotted members of HysteresisModel (non-extendable)
    __slots__ = ("dH", "units", "H_major", "M_major", "remanence", "coercivity")
    
    def __init__(self, dH, units="SI"):
        
        # Raise exceptions for invalid input parameters
        if dH<=0: raise ValueError("\"dH\" (field step-size) must be a positive, non-zero number")
        if units.upper() not in self.UNITS:
            raise ValueError("\"units\" must be one of {:s}".format(', '.join(self.UNITS)))
            
        # Set input parameters
        self.dH = dH
        self.units = self.UNITS[units.upper()]
        
    def _get_major(self, sat_tol):
        """VIRTUAL: computes the major H-M curve & its critical points"""
        raise NotImplemented("This is a virtual method. Stop.")
        
    def _convert_params(self, cf):
        """VIRTUAL: converts model-specific parameter field units"""
        return 0
        
    def path(self, H_path, M0):
        """VIRTUAL: computes hysteresis curves along applied field paths"""
        raise NotImplemented("This is a virtual method. Stop.")
                
    def point(self, H_point, curve="upper"):
        """Computes magnetization on a curve of the major loop at a particular field strength"""
        
        # Select points along the desired curve
        if curve.lower()=='initial': indices = (0, int(.2*len(self.H_major)))
        elif curve.lower()=='upper': indices = (int(.2*len(self.H_major)), int(.6*len(self.H_major)))
        elif curve.lower()=='lower': indices = (int(.6*len(self.H_major)), int(len(self.H_major)))
        else: raise ValueError("\"curve\" must be one of: initial, upper, lower")
        H = self.H_major[indices[0]:indices[1]]
        M = self.M_major[indices[0]:indices[1]]
        
        # Find curve point closest to prescribed point
        H_dists = [abs(h - H_point) for h in H]
        neighbors = argsort(H_dists)[:2]
            
        # Interpolate magnetic flux density between neighboring points
        return sum([M[pt]*H_dists[pt] for pt in neighbors])/self.dH
    
    def convert_units(self, units):
        """Converts model field units"""
        if units.upper() not in self.UNITS:
            raise ValueError("\"units\" must be one of {:s}".format(', '.join(self.UNITS)))
        cf = self.UNITS[units.upper()]/self.units
        self.units *= cf
        self.dH *= cf
        self.H_major = [cf*h for h in self.H_major]
        self.M_major = [cf*m for m in self.M_major]
        self.remanence = [cf*r for r in self.remanence]
        self.coercivity = [cf*c for c in self.coercivity]
        self._convert_params(cf)

    def save(self, path="./ja_model.pkl"):
        """Saves a hysteresis model object as a .pkl file"""
        with open(path,'wb') as file:
            dump(self, file)
    
    @classmethod
    def load(cls, path):
        """Loads a hysteresis model object from a pkl file"""
        with open(path, 'rb') as file:
            return load(file)
        