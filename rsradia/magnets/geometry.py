import numpy as np


def get_circle_center(p1, p2, r):
    """Return the center (x,y) of a circle given two points p1 and p2 and radius r"""
    x1, y1 = p1
    x2, y2 = p2
    
    c = x1**2 + y1**2 - x2**2 - y2**2
    x3 = x1 - x2
    y3 = y1 - y2
    cy = 0.5 * c / y3 - y1
    c2 = (x1**2 + cy**2 - r**2)

    a = (x3**2 + y3**2) / y3**2
    b = -2* x1 - 2 * cy * x3 / y3
    center0 = (-b - np.sqrt(b**2 - 4 * a * c2)) / (2 * a)
    center1 = -x3 * center0 / y3 + 0.5 * c / y3
    
    return center0, center1

def get_arc_points(x1, x2, center, radius, N=15):
    """Return a list of points (x, y) giving an arc on a circle at center=(x, y) and radius. Arc is between points on the x-axis [x1, x2)"""
    x = np.linspace(x1, x2, N)
    
    a = 1.0
    b = -2 * center[1]
    c = center[1]**2  - radius**2 + (x - center[0])**2
    
    y = (- b - np.sqrt(b**2 - 4 * a *c )) / (2 * a)
    
    return x, y 

def get_intersection(line1_p1, line1_p2, line2_p1, line2_p2):
    """find the intersection of two infinite lines defined by two points on each line"""
    x1, y1 = line1_p1
    x2, y2 = line1_p2
    x3, y3 = line2_p1
    x4, y4 = line2_p2
    
    D = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    
    Px = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
    Px /= D
    
    Py = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
    Py /= D
    
    return Px, Py