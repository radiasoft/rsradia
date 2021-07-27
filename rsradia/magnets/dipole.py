import radia as rad


def _create_point_table(pole_width, pole_separation, pole_height, top_height, leg_width, gap_height):
    """
    Construct 2D slice of an H-dipole in the YZ plane.

    All distances are half-widths/lengths of the full magnet. Points are specified in the first quadrant and mirrored
    to other quadrants later. The point list is returned in counter-clockwise order starting from the midpoint above
    the gap center on the pole face.

    :param pole_width: (float)
    :param pole_separation: (float)
    :param pole_height: (float)
    :param top_height: (float)
    :param leg_width: (float)
    :param gap_height: (float)
    :return:
    """
    p1 = [0.0, top_height + pole_height + gap_height]
    p2 = [pole_width + pole_separation, top_height + pole_height + gap_height]

    p3 = [pole_width + pole_separation + leg_width, gap_height + pole_height + top_height]
    p4 = [pole_width + pole_separation + leg_width, 0.0]
    p5 = [pole_width + pole_separation, 0.0]

    p6 = [pole_width + pole_separation, pole_height + gap_height]
    p7 = [pole_width, pole_height + gap_height]
    p8 = [pole_width, gap_height]

    p_middle = [0.0, gap_height]

    point_table = [p1, p2, p3, p4, p5, p6, p7, p8, p_middle]

    return point_table[::-1]


def _get_all_points_top(table):
    """Reflect point table from `create_point_table` to quadrant 2 to form top of H-dipole"""

    coordinates = []
    for point in table:
        if point[0] != 0.:
            reflection = [-point[0], point[1]]
            coordinates.append(reflection)
    coordinates = table + coordinates[::-1]

    return coordinates


def _get_all_points_bottom(table):
    """Reflect point table from `get_all_points_bottom` to quadrant 3 & 4 to form bottom of H-dipole"""
    coordinates = []
    for point in table:
        reflection = [point[0], -point[1]]
        coordinates.append(reflection)

    return coordinates[::-1]


def create_pole(coordinates, center, length, mode=0, triangle_min_size=0.5, triangle_max_size=1.0):
    """
    Form geometry for the full pole piece of an H-dipole using radia.ObjMltExtTri.

    :param coordinates: (list) List of points defining the boundaries of the pole piece in the YZ plane.
    :param center: (float) Center point of dipole in x (longitudinal center for beam frame).
    :param length: (float) Length of the dipole in x
    :param mode: (int) If 0 (default) then the pole piece is divisioned into polygons based on point ordering from
    coordinate list. If != 0 then a Triangular mesh is automatically generated.
    :param triangle_min_size: (float) Only used if mode != 0. Sets the minimum triangle area for automatic division.
    :param triangle_max_size: (float) Only used if mode != 0. Sets the maximum triangle area for automatic division.
    :return: Radia object containing top and bottom pole pieces.
    """
    x = center
    lx = length
    pole_sub = [[1, 1] for _ in range(len(coordinates))]
    # simple build
    if not mode:
        pole = rad.ObjMltExtTri(x, lx, coordinates, pole_sub)
    else:
        str_param = 'ki->Numb,TriAngMin->' + str(triangle_min_size) + ',TriAreaMax->' + str(triangle_max_size)
        pole = rad.ObjMltExtTri(x, lx, coordinates, pole_sub, 'x', [0., 0., 0.], str_param)
    return pole


def make_racetrack_coil(center, radii, sizes, segments=15, current=1):
    """
    Create coil for H-dipole. Coil runs in the XY plane.
    :param center: (list) Center of the coil in [x, y, z].
    :param radii: (list) Inner and outer edges for the coil.
    :param sizes: (list) Straight sections lengths in X and Y; coil height in Z.
    :param segments: (int) Number of segments for coil corners (default: 15).
    :param current: (float) Current carried by the coil (default: 1).
    :return: Radia object representing the coil
    """
    return rad.ObjRaceTrk(center, radii, sizes[:2], sizes[2], segments, current, 'man', 'z')


def make_dipole(pole_dimensions, center, length, current=-10000,
                trimesh_mode=0, longitudinal_divisions=4):
    """
    Construct a complete H-dipole made of iron.
    :param pole_dimensions: (dict) Parameters describing geometry of pole piece. See `_create_point_table`.
    :param center: (float) Center point of dipole in x (longitudinal center for beam frame).
    :param length: (float) Length of the dipole in x
    :param current: (float) Current carried by dipole coils (default: -10000)
    :param trimesh_mode: (int) If 0 (default) then the pole piece is divisioned into polygons based on point ordering
    from coordinate list. If != 0 then a Triangular mesh is automatically generated.
    :param longitudinal_divisions: (int) Number of slices to divide up the dipole into along the x-axis (default: 4)
    :return:
    """
    # coil_factor increases coil size slightly to accommodate sharp corners of pole piece
    coil_length_factor = 1.005
    coil_height_factor = 2. / 3.
    # Geometry for the poles
    table_quadrant_one = _create_point_table(**pole_dimensions)
    top_coodinates = _get_all_points_top(table_quadrant_one)
    bottom_coordinates = _get_all_points_bottom(top_coodinates)

    top_pole = create_pole(top_coodinates, center, length, mode=trimesh_mode)
    bottom_pole = create_pole(bottom_coordinates, center, length, mode=trimesh_mode)

    # Material for the poles (uses Iron)
    ironmat = rad.MatSatIsoFrm([20000, 2], [0.1, 2], [0.1, 2])
    rad.MatApl(top_pole, ironmat)
    rad.MatApl(bottom_pole, ironmat)

    # Coils
    top_coil = make_racetrack_coil(center=[0, 0.0, pole_dimensions['gap_height'] + pole_dimensions['pole_height'] / 2.],
                                   radii=[0.1, 0.4],
                                   sizes=[length * coil_length_factor, 
                                          pole_dimensions['pole_width'] * 2 * coil_length_factor, 
                                          pole_dimensions['pole_height'] * coil_height_factor],
                                   current=current)
    bottom_coil = make_racetrack_coil(center=[0, 0.0, -1. * (pole_dimensions['gap_height'] + pole_dimensions['pole_height'] / 2.)],
                                      radii=[0.1, 0.4],
                                      sizes=[length * coil_length_factor, 
                                             pole_dimensions['pole_width'] * 2 * coil_length_factor, 
                                             pole_dimensions['pole_height'] * coil_height_factor],
                                      current=current)

    # Visualization
    rad.ObjDrwAtr(top_pole, [0, 0.4, 0.8])
    rad.ObjDrwAtr(bottom_pole, [0, 0.4, 0.8])
    rad.ObjDrwAtr(top_coil, [0.2, 0.9, 0.6])
    rad.ObjDrwAtr(bottom_coil, [0.2, 0.9, 0.6])

    # Element Division
    rad.ObjDivMag(top_pole, [longitudinal_divisions, 1, 1])
    rad.ObjDivMag(bottom_pole, [longitudinal_divisions, 1, 1])

    return rad.ObjCnt([top_pole, bottom_pole, top_coil, bottom_coil])