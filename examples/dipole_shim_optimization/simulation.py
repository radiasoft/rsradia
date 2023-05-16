import dipole as dipole_module
import numpy as np
import radia as rad

# Settings for the dipole
DIPOLE_X_CENTER = 0.0
DIPOLE_X_LENGTH = 25.0
DIPOLE_CURRENT = 2000.
BASELINE = {
    'pole_width': 1.0,
    'pole_separation': 1.5,
    'pole_height': 1.5,
    'top_height': 1.5,
    'gap_height': 0.75,
    'leg_width': 1.5
}
OBJ_FILENAME = 'obj_values.npy'


# Modfications to add the shim
def _make_shim_point_table(x_center, y_center, radius, count):
    theta = np.linspace(0, -np.pi, count)
    x, y = radius * np.cos(theta) + x_center, radius * np.sin(theta) + y_center
    
    return np.array([x, y]).T.tolist()
    

def make_dipole_with_shim(radius, dipole_dict, triangle_min_angle, triangle_max_size):
    if radius is not None:
        x_center = dipole_dict['pole_width'] - radius
        shim_table = _make_shim_point_table(x_center, dipole_dict['gap_height'], radius, 25)
    else:
        shim_table = [[dipole_dict['pole_width'], dipole_dict['gap_height']]]
    shim_table.append([0., dipole_dict['gap_height']])
    dipole = dipole_module.make_dipole(dipole_dict, DIPOLE_X_CENTER, DIPOLE_X_LENGTH,
                                       trimesh_mode=1,
                                       triangle_min_size=triangle_min_angle, triangle_max_size=triangle_max_size,
                                       current=DIPOLE_CURRENT, 
                                       pole_face_array=shim_table)
    return dipole


# helper functions used to calculate the value of the objective function
def _bz_on_mesh():
    Y, Z = np.meshgrid(np.linspace(-baseline['pole_width']*0.5, baseline['pole_width']*0.5, 64),
                   np.linspace(-baseline['gap_height']*0.75, baseline['gap_height']*0.75, 64))

    mesh = np.stack([np.zeros_like(Y).ravel(), Y.ravel(), Z.ravel()]).T

    mesh = [list(m) for m in mesh]

    Bz_on_mesh = rad.Fld(dipole, 'Bz', mesh)
    
    return Bz_on_mesh

def _bz_center():
    return rad.Fld(dipole, 'Bz', [0, 0, 0])

def _gfr_on_grid():
    Bz_on_mesh = _bz_on_mesh()
    Bz0 = _bz_center()

    return np.sum(gfr(Bz_on_mesh, Bz0).reshape(*Y.shape) < 1e-4)

def _field_in_plane(y_min, y_max, obj, steps=25):
    int_field = []
    yc = np.linspace(y_min, y_max, steps)
    for y in yc:
        v = rad.FldInt(obj, 'inf', 'ibz', [-150, y, 0.], [150, y, 0.])
        int_field.append(v)
    return np.array([yc, int_field])


def obj_f(obj):
    "Objective function formatted for DFO-LS"
    field_prof = _field_in_plane(-1., 1., obj, steps=125)
    int_bz0 = rad.FldInt(obj, 'inf', 'ibz', [-150, 0., 0.], [150, 0., 0.])
    res = field_prof[1, :] - int_bz0
    
    np.save(OBJ_FILENAME, res)
    
    return np.sum(res**2),  res


def parallel_obj_f(_):
    "Load output of `obj_f` and calculate result. Should be used if simulation is run in parallel (see yaml file for additional information."
    fvec =  np.load(OBJ_FILENAME)
    return np.sum(fvec**2), fvec


def main(radius, precision=1e-4, maxiter=4500, min_angle=20., max_area=0.5):
    dipole = make_dipole_with_shim(radius, BASELINE,
                              triangle_min_angle=min_angle,
                              triangle_max_size=max_area)
    rad.Solve(dipole, precision, maxiter)
    
    np.save('field_in_plane.npy', _field_in_plane(-1., 1., dipole, steps=125))
    
    return obj_f(dipole)
