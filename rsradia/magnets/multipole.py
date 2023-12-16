import radia as rad
import numpy as np
from math import *
from rsradia.magnets import geometry

π = pi

def _pole_tip_path(n_poles, width, gap, height, n_curve = 6,
                  poletip_frac = 0.5, yoketip_frac = 0.6):
    """
    Return a discretized path that outlines the pole tip cross-section.

    Arguments:
        n_poles      = number of magnet pole tips (even integer)
        width        = width of pole tip / mm
        gap          = distance between opposing pole tips / mm
        height       = height of pole tip / mm
        n_curve      = number of intervals for discretizing (half) the pole face
        poletip_frac = lower fraction of pole tip, possibly subject to finer segmentation
        yoketip_frac = ratio of yoke height (or depth) to pole tip width

    The pole tip’s shape is that of a hyperbola:
      z^2 - (hyp * y)^2 = z0^2, which has asymptotes z = ±hyp * y.
    """

    # define some geometric parameters
    tan_np = tan(π / n_poles)
    hyp = 1 / tan_np  # slope of hyperbola's asymptote
    z0 = gap / 2
    ym = width / 2
    zm = hypot(hyp * ym, z0)
    # sanity check: pole tip includes all of pole face?
    assert zm < z0 + height, \
          "Pole tip height too short to accommodate entire curved pole face."

    # construct hyperbolic path
    dy = ym / n_curve
    ny = n_curve + 2  # go two points beyond so we don't have to extend the array
    Γ_tip = [ [iy * dy, hypot(hyp * iy * dy, z0)] for iy in range(ny + 1) ]
    # and
    # modify last two points so as to outline lower portion of the (half) pole tip
    ht_lower = poletip_frac * height
    # sanity check: lower fraction of pole tip includes all of pole face?
    assert zm < z0 + ht_lower, \
          "Lower fraction of pole tip cannot accommodate entire pole face."
    Γ_tip[n_curve + 1] = [ym, z0 + ht_lower]
    Γ_tip[n_curve + 2] = [ 0, z0 + ht_lower]

    return Γ_tip


def create_multipole(n_poles, thick, width, gap, height, chamfer, tip_coil_sep,
                     iron_mat, curr_density,
                     r_min = 2., clearance = 2.,
                     poletip_frac = 0.5, yoketip_frac = 0.6,
                     chamfer_ang = 45., skew = False, n_curve = 6,
                     nx = 2, ny = 2, nzlt = 3, nzut = 3, nat = 2, nycb = 3, nac = 2,
                     iron_color = [0.0, 0.7, 0.9], coil_color = [1.0, 0.3, 0.3]):
    """
    Return a Radia representation of a simple multipole electromagnet with the shape
    of its pole tip modified according to the parameter arrays αl and βl.

    Required Arguments:
        n_poles      = number of magnet pole tips (even integer)
        thick        = length of iron yoke along particle trajectory / mm
        width        = width of pole tip / mm
        gap          = distance between opposing pole tips / mm
        height       = height of pole tip / mm
        chamfer      = size of chamfer on pole tip ends / mm
        tip_coil_sep = distance by which coil edge is set back from the pole tip / mm
        iron_mat     = Radia representation of the iron’s magnetic characteristics
                         (e.g. M-H curve)
        curr_density = current density / (A / mm^2)
    Optional Arguments:
        r_min        = minimum coil radius / mm
        clearance    = clearance between coil corner and diagonal between sectors / mm
        poletip_frac = lower fraction of pole tip, possibly subject to finer segmentation
        yoketip_frac = ratio of yoke height (or depth) to pole tip width
        chamfer_ang  = angle of chamfer normal w/rt pole tip axis / deg
        skew         = False | True | (angle from ‘normal’ / deg)
        n_curve      = number of intervals for discretizing (half) the pole face
        nx           = number of segments along X axis, within distance thick / 2
        ny           = number of segments along Y axis, within distance width / 2
        nzlt         = number of segments along Z axis, lower portion of pole tip
        nzut         = number of segments along Z axis, upper portion of pole tip
        nat          = number of azimuthal segments at top of pole tip
        nycb         = number of segments along Y axis, along the yoke cross bar
        nac          = number of azimuthal segments at corner of yoke
        iron_color   = color to use for iron yoke and pole tips
        coil_color   = color to use for current-carrying coils

    In the above context, the coordinate axes X, Y, Z respectively align with the
    beam trajectory, across the pole tip, and parallel to the pole tip, with the
    origin at the center of the magnet.

    This function constructs one-fourth (right front) of one sector of a multipole
    magnet. It then applies appropriate symmetries to construct the full magnet,
    and then orients the magnet as desired.

    To Check: Does positive current correspond to positive field strength?
              Does skew have the correct orientation?
    """
    # sanity check: even number of magnetic poles?
    assert n_poles % 2 == 0, "Argument n_poles must equal an even integer."
    # sanity check: positive coil height?
    assert tip_coil_sep < height, "Tip-coil separation exceeds height of pole tip."
    # sanity check: chamfers do not cut into all of pole tip?
    assert chamfer < thick / 2, "Chamfer too large."

    # define a few useful vectors
    ctr   = [0, 0, 0]
    x_hat = [1, 0, 0]
    y_hat = [0, 1, 0]
    z_hat = [0, 0, 1]

    # define segmentation parameters
    # :: [nx, ny, nz] or [nr, na, nl]
    n1 = [nx, ny,   nzlt]  # lower pole tip
    n2 = [nx, ny,   nzut]  # upper pole tip
    n3 = [ny, nat,  nx  ]  # top of pole tip
    n4 = [nx, nycb, ny  ]  # cross bar
    n5 = [ny, nac,  nx  ]  # corner

    # define some derived geometric parameters
    tan_np = tan(π / n_poles)
    hyp = 1 / tan_np  # slope of hyperbola's asymptote
    z0 = gap / 2
    ym = width / 2
    zm = hypot(hyp * ym, z0)
    ht_lower = poletip_frac * height

    # create the discretized path that outlines the pole tip cross-section
    Γ_tip = _pole_tip_path(n_poles, width, gap, height, n_curve = 6,
                           poletip_frac = poletip_frac, yoketip_frac = yoketip_frac)

    # create and segment the lower portion of the (half) pole tip
    g_tip = rad.ObjThckPgn(thick / 4, thick / 2, Γ_tip)
    rad.ObjDivMag(g_tip, n1)

    # create and segment the upper portion of the (half) pole tip
    ht_upper = height - ht_lower
    g_pole = rad.ObjRecMag([thick / 4, width / 4, z0 + height - ht_upper / 2],
                       [thick / 2, width / 2, ht_upper])
    rad.ObjDivMag(g_pole, n2)

    # combine the lower and upper portions of the (half) pole tip
    g_pt = rad.ObjCnt([g_tip, g_pole])
    # and
    # cut chamfer, then retain desired metal
    θ = chamfer_ang * degree
    g_poletip = rad.ObjCutMag(g_pt, [thick / 2 - chamfer, 0, z0],
                              [sin(θ), 0, -cos(θ)])[0]

    # create and segment "corner" above (half) pole tip
    depth = yoketip_frac * width
    g_top = rad.ObjRecMag([thick / 4, width / 4, z0 + height + depth / 2],
                       [thick / 2, width / 2, depth])
    cy = [ [[0, ym, z0 + height], x_hat], [0, 0, z0 + height], 2 * depth / width ]
    rad.ObjDivMag(g_top, n3, 'cyl', cy)

    # create and segment horizontal yoke segment to corner
    length = tan_np * (z0 + height) - ym
    g_bar = rad.ObjRecMag([thick / 4, ym + length / 2, z0 + height + depth / 2],
                       [thick / 2,length, depth])
    rad.ObjDivMag(g_bar, n4)

    # outline the corner
    yc = ym + length
    zc = z0 + height
    Γ_corner = [[yc, zc], [yc, zc + depth], [yc + depth * tan_np, zc + depth]]
    # and
    # create and segment yoke corner
    g_corner = rad.ObjThckPgn(thick / 4, thick / 2, Γ_corner)
    cy = [[[0, yc, zc], x_hat], [0, yc, zc + depth], 1]
    rad.ObjDivMag(g_corner, n5, 'cyl', cy)

    # create container for the (half) pole tip plus attached crossbar
    g_yoke = rad.ObjCnt([g_poletip, g_top, g_bar, g_corner])
    # specify the iron
    rad.MatApl(g_yoke, iron_mat)
    # and set color for iron
    rad.ObjDrwAtr(g_yoke, iron_color)

    # create coil1
    ht_coil = height - tip_coil_sep
    # sanity check: coil does not extend below outer edge of curved pole tip
    assert zm < z0 + height - ht_coil, \
           "Inner coil will obscure part of the curved pole tip."
    wd_to_diagonal = (gap / 2 + tip_coil_sep) * tan_np
    r1 = wd_to_diagonal - clearance - ym + r_min
    coil1 = rad.ObjRaceTrk([0, 0, z0 + height - ht_coil / 2], [r_min, r1],
                           [thick, width - 2 * r_min],
                           ht_coil, 3, curr_density)
    # and set color for coil1
    rad.ObjDrwAtr(coil1, coil_color)

    # create coil2
    ht_coil = (height - tip_coil_sep) / 2
    wd_to_diagonal = (z0 + height - ht_coil) * tan_np
    r2 = wd_to_diagonal - clearance - ym + r_min
    coil2 = rad.ObjRaceTrk([0, 0, z0 + height - ht_coil / 2], [r1, r2],
                           [thick, width - 2 * r_min],
                           ht_coil, 3, curr_density)
    # and set color for coil2
    rad.ObjDrwAtr(coil2, coil_color)

    # apply symmetries to create full pole tip plus attached crossbar
    # :: reflection in y-z plane, with zero field perpendicular to the plane
    rad.TrfZerPerp(g_yoke, ctr, x_hat)
    # :: reflection in z-x plane, with zero field perpendicular to the plane
    rad.TrfZerPerp(g_yoke, ctr, y_hat)

    # create container for full magnet: here iron yoke plus coils in one sector
    g_magnet = rad.ObjCnt([g_yoke, coil1, coil2])

    # :: reflection across diagonal plane, with zero field parallel to the plane
    rad.TrfZerPara(g_magnet, ctr, [0, cos(π / n_poles), sin(π / n_poles)])
    # ==>> at this point we have a matched pair of pole tips
    #      they subtend an angle 2 * (2π / n_poles) = 4π / n_poles

    # apply rotation symmetries to create full multipole electromagnet
    rad.TrfMlt(g_magnet, rad.TrfRot(ctr, x_hat, 4 * π / n_poles), int(n_poles / 2))

    # ensure upright orientation of this multipole
    if n_poles % 4 == 0:
        rad.TrfOrnt(g_magnet, rad.TrfRot(ctr, x_hat, π / n_poles))

    # adjust orientation for skew multipole
    if skew == False:
        skew_angle = 0.
    elif skew == True:
        skew_angle = (π / n_poles)
    elif isinstance(skew, numbers.Number):
        skew_angle = skew * degree
    else:
        assert False, "The argument skew must equal one of " \
                      "True, False, or numeric angle in degrees."
    if skew_angle != 0.:
        rad.TrfOrnt(g_magnet, rad.TrfRot(ctr, x_hat, skew_angle))

    return g_magnet
