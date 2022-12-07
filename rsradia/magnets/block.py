### block.py
### November 2022
### A simple rectangular ferromagnet wrapped with a racetrack coil

from numpy import array
from radia import ObjCnt, ObjRecMag, ObjRaceTrk, ObjDrwAtr, Fld, UtiDel
from jupyter_rs_radia import radia_viewer as rv

from .. import MU0
from . import TRI_MIN_SIZE, TRI_MAX_SIZE, UtiDelStuf

class Block():
    '''
    A simple rectangular ferromagnet wrapped with a racetrack coil

    parameters:
        center: cartesian center point of the magnet
        block_dims: dimensions of the central block (l, w, h)
        coil_inner: inner radius of the racetrack curves
        coil_thickness: thickness of the racetrack coil
        coil_gap: gap between racetrack straights & the central block
        coil_divs: coil_divs: number of divisions along coil bent sections
        name: a name to use in RadiaViewer windows
        block_col: color assigned to the central block
        coil_col: color assigned to the racetrack coil
    '''

    def __init__(self, center, block_dims, coil_inner, coil_thickness, coil_gap, coil_divs=15,\
        name="Block Magnet", block_col=[0, 0.4, 0.8], coil_col=[0.2, 0.9, 0.6]):

        # Assign geometric parameters to class members
        self._center = center
        self._block_dims = block_dims
        self._coil_thickness = coil_thickness
        self._coil_inner = coil_inner
        self._coil_gap = coil_gap
        self._coil_divs = coil_divs

        # Assign aesthetic parameters to class members
        self._name = name
        self._block_col = block_col
        self._coil_col = coil_col

        # Initialize a dictionary containing elements of the magnet
        self._build()        

    # Constructs the magnet in a completely unmagnetized state
    def _build(self):

        # Construct a Radia rectangular magnetic object
        block = ObjRecMag(self._center, self._block_dims)

        # Construct a Radia racetrack coil object
        r_coil = [self._coil_inner, self._coil_inner+self._coil_thickness]
        l_coil = [dim+self._coil_gap for dim in self._block_dims[:2]]
        coil = ObjRaceTrk(self._center, r_coil, l_coil, self._block_dims[2], self._coil_divs, 0, 'man', 'z')

        # Assign the block & coil as magnet elements
        self.elements = {
            'block': block,
            'coil': coil
        }

    # Returns the requested Radia-computed field for the magnet (at its center, by default)
    def field(self, typestr, pt=[0,0,0]):
        temp_container = ObjCnt([self.elements['block'], self.elements['coil']])
        out = Fld(temp_container, typestr, pt)
        UtiDel(temp_container)
        return out

    # Magnetizes the magnet by applying a current & accounting for hysteresis
    def magnetize(self, J, hysteresis_model):

        # Store the previous magnetization state (in SI units)
        H_prev = array(Fld(self.elements['coil'], 'h', self._center))/MU0 #array(self.field('h', self._center))/MU0
        M_prev = array(Fld(self.elements['block'], 'm', self._center))/MU0 #array(self.field('m', self._center))/MU0

        # Delete old magnet elements
        UtiDelStuf(self.elements['block'])
        UtiDelStuf(self.elements['coil'])

        # Construct a new coil object with the specified current applied
        r_coil = [self._coil_inner, self._coil_inner+self._coil_thickness]
        l_coil = [dim+self._coil_gap for dim in self._block_dims[:2]]
        self.elements['coil'] = ObjRaceTrk(self._center, r_coil, l_coil, self._block_dims[2], self._coil_divs, J, 'man', 'z')
        H_next = array(Fld(self.elements['coil'], 'h', self._center))/MU0
            
        # Construct a new block object with magnetization specified by the hysteresis model
        M_hyst = [0, 0, 0]
        for i in range(3):
            H_hyst, B_hyst = hysteresis_model.path(array([[H_prev[i], H_next[i]]]), M_prev[i])
            M_hyst[i] = B_hyst[-1]
        self.elements['block'] = ObjRecMag(self._center, self._block_dims, M_hyst)
        
        return H_next, M_hyst

    # Views the magnet using its RadiaViewer member
    def display(self):

        # Set color attributes for magnet elements
        ObjDrwAtr(self.elements['block'], self._block_col)
        ObjDrwAtr(self.elements['coil'], self._coil_col)

        # Create a temporary container for magnet elements & add it to the viewer
        viewer = rv.RadiaViewer()
        temp_container = ObjCnt([self.elements['block'], self.elements['coil']])
        viewer.add_geometry(self._name, temp_container)

        # Open the viewer and delete the temporary container
        display = viewer.display()
        UtiDel(temp_container)
        return display