### block.py
### November 2022
### A simple rectangular ferromagnet wrapped with a racetrack coil

from numpy import zeros
from radia import Fld, ObjCnt, ObjRecMag, ObjRaceTrk, ObjDrwAtr, UtiDel
from jupyter_rs_radia import radia_viewer as rv

from .. import PI, MU0
from . import TRI_MIN_SIZE, TRI_MAX_SIZE

class Block():
    '''
    A simple rectangular ferromagnet wrapped with a racetrack coil

    parameters:
        center: cartesian center point of the magnet
        block_dims: dimensions of the central block (l, w, h)
        coil_thickness: thickness of the racetrack coil
        coil_inner: inner radius of the racetrack curves
        coil_gap: gap between racetrack straights & the central block
        name: a name to use in RadiaViewer windows
        block_col: color assigned to the central block
        coil_col: color assigned to the racetrack coil
    '''

    def __init__(self, center, block_dims, coil_inner, coil_thickness, coil_gap, \
        name="Block Magnet", block_col=[0, 0.4, 0.8], coil_col=[0.2, 0.9, 0.6]):

        # Assign geometric parameters to class members
        self._center = center
        self._block_dims = block_dims
        self._coil_thickness = coil_thickness
        self._coil_inner = coil_inner
        self._coil_gap = coil_gap

        # Assign aesthetic parameters to class members
        self._name = name
        self._block_col = block_col
        self._coil_col = coil_col
        self._viewer = rv.RadiaViewer()

        # Initialize a dictionary containing elements of the magnet
        self._build()        

    # Constructs the magnet in a completely unmagnetized state
    def _build(self):

        # Construct a Radia rectangular magnetic object
        block = ObjRecMag(self._center, self._cube_dims)

        # Construct a Radia racetrack coil object
        r_coil = [self._coil_inner, self._coil_inner+self._coil_thickness]
        coil = ObjRaceTrk(self._center, r_coil, self._block_dims[:2], self._block_dims[2], 15, 0, 'man', 'z')

        # Assign the block & coil as magnet elements
        self._elements = {
            'block': block,
            'coil': coil
        }

    # Returns the requested Radia-computed field for the magnet (at its center, by default)
    def field(self, typestr, pt=self._center):
        temp_container = ObjCnt([self._elements['block'], self._elements['coil']])
        out = Fld(temp_container, typestr, pt)
        UtiDel(temp_container)
        return out

    # Magnetizes the magnet by applying a current & accounting for hysteresis
    def magnetize(self, J, hysteresis_model):

        # Store the previous magnetization state (in SI units)
        H_prev = self.field('h', self._center)/MU0
        M_prev = self.field('m', self._center)/MU0

        # Delete old magnet elements
        UtiDel(self._elements['block'])
        UtiDel(self._elements['coil'])

        # Construct a new coil object with the specified current applied
        r_coil = [self._coil_inner, self._coil_inner+self._coil_thickness]
        self._elements['coil'] = ObjRaceTrk(self._center, r_coil, self._block_dims[:2], self._block_dims[2], 15, J, 'man', 'z')
        H_next = array(Fld(coil, 'h', self._center))/MU0
            
        # Construct a new block object with magnetization specified by the hysteresis model
        Mhyst = [0, 0, 0]
        for i in range(3):
            Hhyst, Bhyst = hysteresis_model.path(array([[H_prev[i], H_next[i]]]), M_prev)
            Mhyst[i] = Bhyst[-1]
        self._elements['block'] = ObjRecMag(self._center, self._cube_dims, Mhyst)

    # Views the magnet using its RadiaViewer member
    def view(self):

        # Set color attributes for magnet elements
        ObjDrwAtr(self._elements['block'], self._block_col)
        ObjDrwAtr(self._elements['coil'], self._coil_col)

        # Create a temporary container for magnet elements & add it to the viewer
        temp_container = ObjCnt([self._elements['block'], self._elements['coil']])
        _viewer.add_geometry(self._name, self._container)

        # Open the viewer and delete the temporary container
        _viewer.display()
        UtiDel(temp_container)