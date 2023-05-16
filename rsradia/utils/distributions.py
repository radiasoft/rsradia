### distributions.py
### July 2022
### Distribution functions for probabilities defined over 2D grids

from numpy import array, zeros, exp, pi, einsum
from scipy.linalg import det, inv

# Note: in all functions below, the "grid" argument should be a numpy array with
# shape (n, 2), in which n is the total number of pixels (n = nx * ny) and the 
# 2 elements per pixel are its x/y values in the grid.

def UNIFORM(grid):
    """A discrete uniform grid density"""
    n_grid = array(grid.shape).prod()
    return zeros(n_grid)+1/n_grid

def GAUSSIAN(grid, mu=array([0, 0]), Sigma=array([[.1, 0], [0, .1]])):
    """A discrete Gaussian density"""
    
    # ARGS
    # mu: distribution mean (x/y coordinate pair)
    # Sigma: distribution covariance matrix
    
    # Determine distances between grid points & the mean
    dists = grid-mu

    # Looping over mixands, add to mixture pdf
    reg = ((2*pi)**2 * det(Sigma))**-.5
    pX = reg*exp(-0.5*einsum('ij,jl,il->i', dists, inv(Sigma), dists))

    return pX/pX.sum()

def GMIXTURE(grid, mus=array([[0, 0]]), Sigmas=array([[[.1, 0], [0, .1]]]), ws=[1]):
    """A distrete Gaussian mixture density"""
    
    # ARGS
    # mus: distribution means (x/y coordinate pairs) for each mixand
    # Sigmas: distribution covariance matrices for each mixand
    # ws: weights for each mixand

    # Define number of mixands & dimensions
    N, k = mus.shape

    # Initialize the mixture pdf
    pX = zeros(grid.shape[0])

    # Looping over mixands, add to mixture pdf
    for n in range(N):
        reg = ((2*pi)**k*det(Sigmas[n]))**-.5
        dists = (grid-mus[n])
        pX += ws[n]*reg*exp(-0.5*einsum('ij,jl,il->i', dists, inv(Sigmas[n]), dists))

    return pX/pX.sum()