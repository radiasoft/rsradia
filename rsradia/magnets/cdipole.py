
import radia as rad

def make_cdipole(gap, thick, width, chamfer, lpp, lap, lha, car, iron,
                  rmin, rmax, current,
                  n1, n2, n3, n4, n5, n6, ncr, nca, nsc,
                  circ = True):
    """
    create a simple dipole electromagnet
    arguments:
        gap     = distance between iron pole tips / mm
        thick   = thickness of iron pole tip (along particle trajectory) / mm
        width   = width of iron pole tip (transverse to particle trajectory) / mm
        chamfer = size of chamfer on iron pole tip / mm
        lpp     = length of pole piece / mm
        lap     = length of segment above the pole piece / mm
        lha     = length of horizontal arm (between the corners) / mm
        car     = corner aspect ratio
        iron    = material specification of the iron
        rmin    = minimum coil radius / mm
        rmax    = maximum coil radius / mm
        current = coil current / A
        n1      = segmentation parameters for pole tip
        n2      = segmentation parameters for vertical segment above pole tip
        n3      = segmentation parameters for corner above pole tip
        n4      = segmentation parameters for horizontal segment between corners
        n5      = segmentation parameters for other corner
        n6      = segmentation parameters for vertical segment inside coil
        ncr     = number of radial segments in the iron corners
        nca     = number of azimuthal segments in the iron corners
        nsc     = number segments in coil corners
        circ    = whether to use cicular segmentation at iron corners (boolean)
    return Radia representation of a simple dipole steering magnet
    """
    # global quantities
    global g1, g2, g3, g4, g5, g6
    global yoke, coil

    # debug parameter
    eps = 0

    # colors
    iron_color = [0.0, 1.0, 1.0]
    coil_color = [1.0, 0.0, 0.0]

    # yoke geometry

    # -- pole piece
    lx1 = thick / 2
    ly1 = width
    lz1 = lpp  # height of pole piece above pole tip / mm
    l1  = [lx1, ly1, lz1]
    k1  = [[thick / 4. - chamfer / 2., 0, gap / 2.],
           [thick / 2. - chamfer, ly1 - 2. * chamfer]]
    k2  = [[thick / 4., 0., gap / 2. + chamfer], [thick / 2., ly1]]
    k3  = [[thick / 4., 0., gap / 2. + lz1],     [thick / 2., ly1]]
    g1  = rad.ObjMltExtRtg([k1, k2, k3])
    rad.ObjDivMag(g1, n1)

    # -- vertical segment above pole piece
    lx2 = thick / 2
    ly2 = ly1
    lz2 = lap
    l2  = [lx2, ly2, lz2]
    p2  = [thick / 4, 0, gap / 2 + lz1 + lz2 / 2 + eps]
    g2  = rad.ObjRecMag(p2, l2)
    rad.ObjDivMag(g2, n2)

    # -- corner above pole piece
    lx3 = thick / 2
    ly3 = ly2
    lz3 = car * ly2
    l3  = [lx3, ly3, lz3]
    p3  = [thick / 4, 0, gap / 2 + lz1 + lz2 + lz3 / 2 + 2 * eps]
    g3  = rad.ObjRecMag(p3, l3)
    #    parameters for circular segmentation
    typ = [[p3[0], p3[1] + ly3 / 2, p3[2] - lz3 / 2], [1, 0, 0],
           [p3[0], p3[1] - ly3 / 2, p3[2] - lz3 / 2], lz3 / ly3]
    #   circular or rectangular segmentation
    if circ:
        rad.ObjDivMag(g3, [ncr, nca, n3[0]], 'cyl', typ)
    else:
        rad.ObjDivMag(g3, n3)

    # -- horizontal segment between corners
    lx4 = thick / 2
    ly4 = lha
    lz4 = lz3
    l4  = [lx4, ly4, lz4]
    p4  = [thick / 4, ly3 / 2 + ly4 / 2 + eps, p3[2]]
    g4  = rad.ObjRecMag(p4, l4)
    rad.ObjDivMag(g4, n4)

    # -- other corner (above coil)
    lx5 = thick / 2
    ly5 = car * lz4
    lz5 = lz4
    l5  = [lx5, ly5, lz5]
    p5  = [thick / 4, p4[1] + (ly4 + ly5) / 2 + eps, p4[2]]
    g5  = rad.ObjRecMag(p5, l5)
    #    parameters for circular segmentation
    typ = [[p5[0], p5[1] - ly5 / 2, p5[2] - lz5 / 2], [1, 0, 0],
           [p5[0], p5[1] + ly5 / 2, p5[2] - lz5 / 2], lz5 / ly5]
    #    circular or square segmentation
    if circ:
        rad.ObjDivMag(g5, [ncr, nca, n5[0]], 'cyl', typ)
    else:
        rad.ObjDivMag(g5, n5)

    # -- vertical segment inside coil
    lx6 = thick/2
    ly6 = ly5
    lz6 = gap/2 + lz1 + lz2
    l6  = [lx6, ly6, lz6]
    p6  = [thick/4, p5[1], p5[2] - (lz6 + lz5)/2 - eps]
    g6  = rad.ObjRecMag(p6, l6)
    rad.ObjDivMag(g6, n6)

    # group iron pieces into a single yoke, set material properties, ...
    yoke = rad.ObjCnt([g1, g2, g3, g4, g5, g6])
    rad.MatApl(yoke, iron)
    # and set color
    rad.ObjDrwAtr(yoke, iron_color)

    # apply symmetry conditions
    rad.TrfZerPerp(yoke, [0, 0, 0], [1, 0, 0])  # across y-z plane, with B parallel to plane
    rad.TrfZerPara(yoke, [0, 0, 0], [0, 0, 1])  # across x-y plane, with B perpendicular to plane

    # coil geometry, current, ...
    hc = 2 * lz6 - rmin
    curr_dens = current / (hc * (rmax - rmin))
    pc = [0, p6[1], 0]
    coil = rad.ObjRaceTrk(pc, [rmin, rmax], [thick, ly6], hc, nsc, curr_dens)
    # and color
    rad.ObjDrwAtr(coil, coil_color)

    # group yoke and coil together, and return
    return rad.ObjCnt([yoke, coil])
