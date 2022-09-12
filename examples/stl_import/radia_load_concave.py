# Example of what happens when you try to load an STL file with concave polyhedron included into Radia

import radia
import os
import trimesh
import time

DIR = '.'
FILE = 'shapr3d_export_2021-09-08_13h01m.stl'
mesh = trimesh.load(os.path.join(DIR, FILE))


radia_object = []

reindex_faces = mesh.faces + 1
vertices_resize = mesh.vertices * 1
print('Starting mesh load')
time.sleep(0.01)
new_obj = radia.ObjPolyhdr(vertices_resize.tolist(), reindex_faces.tolist())
radia_object.append(new_obj)

print("Completed Loading to Radia")