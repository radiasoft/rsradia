import numpy as np
import plotly.graph_objects as go
import trimesh

def plot_mesh(mesh, wireframe=True, incenters=False, nonconvex_faces=False):
    """Plot a mesh from Trimesh"""
    data = []
    
    if wireframe is not None:
        data += [d for d in generate_scatter_data(mesh)]
    if not wireframe:
        data += [d for d in generate_mesh_data(mesh)]
        
    if incenters == 'point':
        data += [d for d in generate_incenters(mesh)]
    elif incenters == 'normal':
        data += [d for d in generate_face_normals(mesh)]
    else:
        pass
    
    if nonconvex_faces:
        data += [d for d in generate_nonconvex_vertices(mesh)]
    
    fig = go.Figure(data=data)
    
    fig.update_layout(
        autosize=False,
        width=1000,
        height=900,
        margin=dict(
            l=50,
            r=50,
            b=100,
            t=100,
            pad=4
        ),

    )
                          
    fig.show()


def extract_lines_from_faces(x, y, z, i, j, k):
    triangles = np.vstack((i,j,k)).T
    vertices = np.vstack((x,y,z)).T
    tri_points = vertices[triangles]

    #extract the lists of x, y, z coordinates of the triangle vertices and connect them by a line
    Xe = []
    Ye = []
    Ze = []
    for T in tri_points:
        Xe.extend([T[k % 3][0] for k in range(4)] + [ None])
        Ye.extend([T[k % 3][1] for k in range(4)] + [ None])
        Ze.extend([T[k % 3][2] for k in range(4)] + [ None])
        
    return Xe, Ye, Ze


def generate_scatter_data(mesh, color='rgb(70,70,70)'):
    if isinstance(mesh, trimesh.Trimesh):
        mesh = [mesh, ]
    for m in mesh:
        x, y, z = m.vertices.T
        i, j, k = m.faces.T
        Xe, Ye, Ze = extract_lines_from_faces(x, y , z, i, j, k)
        yield go.Scatter3d(
                     x=Xe,
                     y=Ye,
                     z=Ze,
                     mode='lines',
                     showlegend=False,
                     line=dict(color=color, width=2)
        )


def generate_mesh_data(mesh, color='rgb(220,220,220)'):
    if isinstance(mesh, trimesh.Trimesh):
        mesh = [mesh, ]
    for m in mesh:
        x, y, z = m.vertices.T
        i, j, k = m.faces.T
        yield     go.Mesh3d(
                        x=x,
                        y=y,
                        z=z,
                        color=color,
                        lighting=dict(ambient=0.3, diffuse=0.75, roughness = 0.4, specular=0.00, fresnel=0.0),
                        i=i,
                        j=j,
                        k=k
        )

def generate_nonconvex_vertices(mesh, color='rgb(88,10,10)'):
    if isinstance(mesh, trimesh.Trimesh):
        mesh = [mesh, ]
    for m in mesh:
        vertex_indices = m.face_adjacency_edges[~m.face_adjacency_convex]
        for edge in vertex_indices:
            Xe, Ye, Ze = [], [], []
            Xe.extend(m.vertices[edge, 0])
            Ye.extend(m.vertices[edge, 1])
            Ze.extend(m.vertices[edge, 2])
            
            yield go.Scatter3d(
                         x=Xe,
                         y=Ye,
                         z=Ze,
                         mode='lines',
                         showlegend=False,
                         line=dict(color=color, width=6)
            )

        
def get_triangle_incenter(vertex1, vertex2, vertex3):
    A, B, C = vertex1, vertex2, vertex3
    a = np.linalg.norm(B - C)
    b = np.linalg.norm(C - A)
    c = np.linalg.norm(A - B)
    
    n = a + b + c
    
    x = a * A / n
    y = b * B / n
    z = c * C / n
    
    return x + y + z
    
        
def generate_incenters(mesh):
    if isinstance(mesh, trimesh.Trimesh):
        mesh = [mesh, ]
        
    X, Y, Z = [], [], []
    for m in mesh:
        for face in m.faces:
            n1, n2, n3 = get_triangle_incenter(*m.vertices[face, :])
            X.append(n1)
            Y.append(n2)
            Z.append(n3)
        yield go.Scatter3d(
                     x=X,
                     y=Y,
                     z=Z,
                     mode='markers',
                     line=dict(color= 'rgb(10,70,30)', width=2)
        ) 
        
def generate_face_normals(mesh):
    if isinstance(mesh, trimesh.Trimesh):
        mesh = [mesh, ]
    color1 = 'blue'
    color2 = 'red'
    X, Y, Z = [], [], []
    
    for m in mesh:
        for face, normal in zip(m.faces, m.face_normals):
            n1, n2, n3 = get_triangle_incenter(*m.vertices[face, :])
            X.append(n1)
            Y.append(n2)
            Z.append(n3)
        yield go.Cone(
                     x=X,
                     y=Y,
                     z=Z,
                     u=m.face_normals[:, 0],
                     v=m.face_normals[:, 1],
                     w=m.face_normals[:, 2],
                     colorscale=[[0, 'rgb(10,70,30)'], [1, 'rgb(10,70,30)']],
                     showscale=False
        ) 

def get_triangle_incenter(vertex1, vertex2, vertex3):
    A, B, C = vertex1, vertex2, vertex3
    a = np.linalg.norm(B - C)
    b = np.linalg.norm(C - A)
    c = np.linalg.norm(A - B)
    
    n = a + b + c
    
    x = a * A / n
    y = b * B / n
    z = c * C / n
    
    return x + y + z
        