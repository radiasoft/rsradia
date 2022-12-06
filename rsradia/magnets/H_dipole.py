### H_dipole.py
### November 2022
### An H-dipole with a ferromagnetic yoke wrapped with racetrack coils

from numpy import array, zeros
from radia import Fld, ObjCnt, ObjMltExtTri, ObjDivMag, ObjRaceTrk, ObjDrwAtr, ObjM, ObjSetM, ObjGeoVol, UtiDel
from jupyter_rs_radia import radia_viewer as rv

from .. import MU0
from . import TRI_MIN_SIZE, TRI_MAX_SIZE, ObjCntStufRec, UtiDelStuf

class HDipole:
    '''
    An H-dipole steering magnet with a ferromagnetic yoke wrapped with racetrack coils

    parameters:
        x_center: center point of the dipole in the x-direction
        length: length of the dipole in the x-direction
        pole_height: height of one-half of the yoke
        pole_width: width of one-half of the yokke
        pole_sep: horizontal separation between poles and yoke body
        gap_height: vertical separation between poles and magnet center
        top_height: height of yoke above/below poles
        leg_width: width
        coil_inner: the inner radius of each racetrack coil
        coil_lfactor: coil length relative to yoke dimensions
        coil_hfactor: coil height relative to yoke dimensions
        coil_rfactor: coil outer radius relative to yoke dimensions
        bevel_base: height of bevel on the yoke base
        x_div: Number of divisions to make in the x-direction
        mesh_mode: 0 -> polygon-based mesh, 1 -> triangle-based mesh
        tri_min_size: minimum allowable size for triangular mesh elements
        tri_max_size: maimum allowable size for triangular mesh elements
        yoke_col: color assigned to the yokes
        coil_col: color assigned to the racetrack coils
    '''

    def __init__(self, x_center, length, pole_height, pole_width, pole_sep, gap_height, top_height, leg_width, coil_inner,
                 coil_lfactor=1.005, coil_hfactor=.75, coil_rfactor=0.85, bevel_base=0.0, x_div=4, coil_div=15, mesh_mode=0,
                 name='H-Dipole', yoke_col=[0, 0.4, 0.8], coil_col=[0.2, 0.9, 0.6]):
        
        # Assign geometric parameters to class members
        self._x_center = x_center
        self._length = length
        self._pole_height = pole_height
        self._pole_width = pole_width
        self._pole_sep = pole_sep
        self._gap_height = gap_height
        self._coil_inner = coil_inner
        self._coil_lfactor = coil_lfactor
        self._coil_hfactor = coil_hfactor
        self._coil_rfactor = coil_rfactor
        self._x_div = x_div
        self._coil_div = coil_div
        self._mesh_mode = mesh_mode
        
        # Define the point table for the magnet
        q1_pts = array([[0.0, gap_height],
                 [pole_width - bevel_base, gap_height],
                 [pole_width, pole_height + gap_height],
                 [pole_width + pole_sep, pole_height + gap_height],
                 [pole_width + pole_sep, 0.0],
                 [pole_width + pole_sep + leg_width, 0.0],
                 [pole_width + pole_sep + leg_width, gap_height + pole_height + top_height],
                 [pole_width + pole_sep, top_height + pole_height + gap_height],
                 [0.0, top_height + pole_height + gap_height]])        
        q2_pts = q1_pts[1:-1]*array([-1, 1])
        self._top_pts = q1_pts.tolist() + q2_pts[::-1].tolist()
        self._bottom_pts = (array(self._top_pts)*array([1,-1])).tolist()      
                
        # Assign aesthetic parameters to class members
        self._name = name
        self._yoke_col = yoke_col
        self._coil_col = coil_col
                
        # Initialize a dictionary containing elements of the magnet
        self._build()
        
    # Constructs the magnet in a completely unmagnetized state
    def _build(self):
        
        # Construct the yoke & divide it into a mesh
        pole_sub = [[1, 1] for _ in range(len(self._top_pts))]
        if not self._mesh_mode:
            top_pole = ObjMltExtTri(self._x_center, self._length, self._top_pts, pole_sub)
            bottom_pole = ObjMltExtTri(self._x_center, self._length, self._bottom_pts, pole_sub)
        else:
            str_param = 'ki->Numb,TriAngMin->' + str(TRI_MIN_SIZE) + ',TriAreaMax->' + str(TRI_MAX_SIZE)
            top_pole = ObjMltExtTri(self._x_center, self._length, self._top_pts, pole_sub, 'x', [0.,0.,0.], str_param)
            bottom_pole = ObjMltExtTri(self._x_center, self._length, self._bottom_pts, pole_sub, 'x', [0.,0.,0.], str_param)
        ObjDivMag(top_pole, [self._x_div, 1, 1])
        ObjDivMag(bottom_pole, [self._x_div, 1, 1])
        
        # Get the volumes of each constituent yoke element
        tp_items = ObjCntStufRec(top_pole)
        bp_items = ObjCntStufRec(bottom_pole)
        self._sub_volumes = zeros((len(tp_items)+len(bp_items)))
        for i in range(len(tp_items)):
            self._sub_volumes[i] = ObjGeoVol(tp_items[i])
        for i in range(len(bp_items)):
            self._sub_volumes[i+len(tp_items)] = ObjGeoVol(bp_items[i])
        
        # Construct the coils
        coil_center = array([0., 0., self._gap_height+self._pole_height/2.])
        radii = array([self._coil_inner, self._pole_sep*self._coil_rfactor])
        sizes = array([self._length*self._coil_lfactor, self._pole_width*2*self._coil_lfactor, self._pole_height*self._coil_hfactor])
        top_coil = ObjRaceTrk(list(coil_center), radii, sizes[:2], sizes[2], self._coil_div, 0, 'man', 'z')
        bottom_coil = ObjRaceTrk(list(-coil_center), radii, sizes[:2], sizes[2], self._coil_div, 0, 'man', 'z')
                    
        # Assign poles & coils as magnet elements
        self._elements = {
            'top_pole': top_pole,
            'top_coil': top_coil,
            'bottom_pole': bottom_pole,
            'bottom_coil': bottom_coil,
        }
                        
    # Returns the requested Radia-computed field for the magnet (at its center, by default)
    def field(self, typestr, pt=[0,0,0]):
        temp_container = ObjCnt([el for el in self._elements.values()])
        out = Fld(temp_container, typestr, pt)
        UtiDel(temp_container)
        return out

    # Magnetizes the magnet by applying a current & accounting for hysteresis
    def magnetize(self, J, hysteresis_model):
                
        # Store the previous magnetization state (in SI units)
        tp_items = ObjCntStufRec(self._elements['top_pole'])
        bp_items = ObjCntStufRec(self._elements['bottom_pole'])
        H_prev = zeros((len(tp_items)+len(bp_items), 3))
        M_prev = zeros((len(tp_items)+len(bp_items), 3))
        for i in range(len(tp_items)):
            vol_center, M = ObjM(tp_items[i])
            H_tc = Fld(self._elements['top_coil'], 'h', vol_center)
            H_bc = Fld(self._elements['bottom_coil'], 'h', vol_center)
            H_prev[i] = array(H_tc)+array(H_bc)
            M_prev[i] = M
        for i in range(len(bp_items)):
            vol_center, M = ObjM(bp_items[i])
            H_tc = Fld(self._elements['top_coil'], 'h', vol_center)
            H_bc = Fld(self._elements['bottom_coil'], 'h', vol_center)
            H_prev[i+len(tp_items)] = array(H_tc)+array(H_bc)
            M_prev[i+len(tp_items)] = M
        H_prev /= MU0
        M_prev /= MU0

        # Delete old magnet elements
        for el in self._elements.values(): UtiDelStuf(el)
        
        # Construct new coil objects with the specified current applied
        coil_center = array([0., 0., self._gap_height+self._pole_height/2.])
        radii = array([self._coil_inner, self._pole_sep*self._coil_rfactor])
        sizes = array([self._length*self._coil_lfactor, self._pole_width*2*self._coil_lfactor, self._pole_height*self._coil_hfactor])
        top_coil = ObjRaceTrk(list(coil_center), radii, sizes[:2], sizes[2], self._coil_div, J, 'man', 'z')
        bottom_coil = ObjRaceTrk(list(-coil_center), radii, sizes[:2], sizes[2], self._coil_div, J, 'man', 'z')
        
        # Construct a new yoke & divide it into a mesh
        pole_sub = [[1, 1] for _ in range(len(self._top_pts))]
        if not self._mesh_mode:
            top_pole = ObjMltExtTri(self._x_center, self._length, self._top_pts, pole_sub)
            bottom_pole = ObjMltExtTri(self._x_center, self._length, self._bottom_pts, pole_sub)
        else:
            str_param = 'ki->Numb,TriAngMin->' + str(TRI_MIN_SIZE) + ',TriAreaMax->' + str(TRI_MAX_SIZE)
            top_pole = ObjMltExtTri(self._x_center, self._length, self._top_pts, pole_sub, 'x', [0.,0.,0.], str_param)
            bottom_pole = ObjMltExtTri(self._x_center, self._length, self._bottom_pts, pole_sub, 'x', [0.,0.,0.], str_param)
        ObjDivMag(top_pole, [self._x_div, 1, 1])
        ObjDivMag(bottom_pole, [self._x_div, 1, 1])
        
        # Assign correct magnetizations to each constituent yoke element
        tp_items = ObjCntStufRec(top_pole)
        bp_items = ObjCntStufRec(bottom_pole)
        M_vol = zeros((len(tp_items)+len(bp_items), 3))
        for i in range(len(tp_items)):
            vol_center, _ = ObjM(tp_items[i])
            H_tc = Fld(top_coil, 'h', vol_center)
            H_bc = Fld(bottom_coil, 'h', vol_center)
            H_vol = array(H_tc)+array(H_bc)
            M_hyst = [0, 0, 0]
            for j in range(3):
                H_hyst, B_hyst = hysteresis_model.path(array([[H_prev[i,j], H_vol[j]]]), M_prev[i,j])
                M_hyst[j] = B_hyst[-1]
            ObjSetM(tp_items[i], M_hyst)
            M_vol[i] = M_hyst
        for i in range(len(bp_items)):
            vol_center, _ = ObjM(bp_items[i])
            H_tc = Fld(top_coil, 'h', vol_center)
            H_bc = Fld(bottom_coil, 'h', vol_center)
            H_vol = array(H_tc)+array(H_bc)
            M_hyst = [0, 0, 0]
            for j in range(3):
                H_hyst, B_hyst = hysteresis_model.path(array([[H_prev[i+len(tp_items),j], H_vol[j]]]), M_prev[i+len(tp_items),j])
                M_hyst[j] = B_hyst[-1]
            ObjSetM(bp_items[i], M_hyst)
            M_vol[i+len(tp_items)] = M_hyst
            
        # Assign poles & coils as magnet elements
        self._elements = {
            'top_pole': top_pole,
            'top_coil': top_coil,
            'bottom_pole': bottom_pole,
            'bottom_coil': bottom_coil,
        }
            
        # Compute the center-point applied field & volume-averaged magnetization
        H_tc = Fld(self._elements['top_coil'], 'h', [self._x_center, 0, 0])
        H_bc = Fld(self._elements['bottom_coil'], 'h', [self._x_center, 0, 0])
        H_next = array(H_tc)+array(H_bc)
        M_next = (self._sub_volumes@M_vol)/self._sub_volumes.sum()
        
        return H_next, M_next

    # Views the magnet using its RadiaViewer member
    def display(self):

        # Set color attributes for magnet elements
        ObjDrwAtr(self._elements['top_pole'], self._yoke_col)
        ObjDrwAtr(self._elements['top_coil'], self._coil_col)
        ObjDrwAtr(self._elements['bottom_pole'], self._yoke_col)
        ObjDrwAtr(self._elements['bottom_coil'], self._coil_col)

        # Create a temporary container for magnet elements & add it to the viewer
        viewer = rv.RadiaViewer()
        temp_container = ObjCnt([el for el in self._elements.values()])
        viewer.add_geometry(self._name, temp_container)

        # Open the viewer and delete the temporary container
        display = viewer.display()
        UtiDel(temp_container)
        return display