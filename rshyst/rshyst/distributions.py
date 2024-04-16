"""distributions.py

A suit of distribution functions used by the rshyst package

LICENSE STATEMENT

Copyright 2024 RadiaSoft LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np
from scipy.linalg import det, inv

# Note: in all functions below, the "grid" argument should be a numpy array with
# shape (n, 2), in which n is the total number of pixels (n = nx * ny) and the 
# 2 elements per pixel are its x/y values in the grid.

def UNIFORM(grid):
    """A discrete uniform grid density"""

    n_grid = np.array(grid.shape).prod()
    
    return np.zeros(n_grid)+1/n_grid

def GAUSSIAN(grid, mu=np.array([0, 0]), Sigma=np.array([[.1, 0], [0, .1]])):
    """A discrete Gaussian density
    
    Args:
      * mu: distribution mean (x/y coordinate pair)
      * Sigma: distribution covariance matrix
    """
        
    # Determine distances between grid points & the mean
    dists = grid-mu

    # Looping over mixands, add to mixture pdf
    reg = ((2*np.pi)**2 * det(Sigma))**-.5
    pX = reg*np.exp(-0.5*np.einsum('ij,jl,il->i', dists, inv(Sigma), dists))

    return pX/pX.sum()

def GMIXTURE(grid, mus=np.array([[0, 0]]), Sigmas=np.array([[[.1, 0], [0, .1]]]), ws=[1]):
    """A distrete Gaussian mixture density
    
    Args:
      * mus: distribution means (x/y coordinate pairs) for each mixand
      * Sigmas: distribution covariance matrices for each mixand
      * ws: weights for each mixand
    """

    # Define number of mixands & dimensions
    N, k = mus.shape

    # Initialize the mixture pdf
    pX = np.zeros(grid.shape[0])

    # Looping over mixands, add to mixture pdf
    for n in range(N):
        reg = ((2*np.pi)**k*det(Sigmas[n]))**-.5
        dists = (grid-mus[n])
        pX += ws[n]*reg*np.exp(-0.5*np.einsum('ij,jl,il->i', dists, inv(Sigmas[n]), dists))

    return pX / pX.sum()
