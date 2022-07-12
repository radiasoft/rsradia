### density_functions.py
### Morgan Henderson, July 2022
### A collection of probability density functions

from numpy import zeros, exp, pi, einsum
from scipy.linalg import det, inv

# Defines a discrete Gaussian pdf p(r) given a mean vector & covariance matrix
def gaussian(r,mu,Sigma):

    # Define number of mixands & dimensions
    k = len(mu)

    # Looping over mixands, add to mixture pdf
    reg = ((2*pi)**k*det(Sigma))**-.5
    pX = reg*exp(-0.5*einsum('ij,jl,il->i',r-mu,inv(Sigma),r-mu))

    return pX/pX.sum()

# Defines a discrete Gaussian mixture pdf p(r) given means, covariances, & weights
def gaussian_mixture(r,mus,Sigmas,ws):

    # Define number of mixands & dimensions
    N,k = mus.shape

    # Initialize the mixture pdf
    pX = zeros(r.shape[0])

    # Looping over mixands, add to mixture pdf
    for n in range(N):
        reg = ((2*pi)**k*det(Sigmas[n]))**-.5
        diff = (r-mus[n])
        pX += ws[n]*reg*exp(-0.5*einsum('ij,jl,il->i',diff,inv(Sigmas[n]),diff))

    return pX/pX.sum()