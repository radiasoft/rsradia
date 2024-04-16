
import radia as rad
from rsradia.magnets import base

def make_block(center, core_dims, coil_inner, coil_thickness, coil_gap, radmat=base.IRONMAT, coil_divs=15):
    
    core = rad.ObjRecMag(center, core_dims)
    r_coil = [coil_inner, coil_inner + coil_thickness]
    l_coil = [dim + coil_gap for dim in core_dims[:2]]
    coil = rad.ObjRaceTrk(center, r_coil, l_coil, core_dims[2], coil_divs, 0, 'man', 'z')

    rad.MatApl(core, radmat)

    return rad.ObjCnt([core, coil])
