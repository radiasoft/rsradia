"""integrators.py

A suit of numerical integrators used by the rshyst package

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

def _get_size(bounds, dt):
    """Determines size of integrator outputs

    Args:
      * bounds: Bounds of integration
      * dt: Step-size over integration domain

    Returns:
      * N: Number of points along integration domain
    """

    N = int(round(np.ceil(abs(bounds[1] - bounds[0])) / abs(dt))) + 1

    return N

def KRON15(fun, bounds, nsplit=1, fun_args={}):
    """A 15-point Gauss-Kronrod quadrature integrator

    Args:
      * fun: Function to be integrated
      * bounds: Bounds of integration
      * nsplit: Number of splits in integration domain (default 1, no splitting)
    """
    
    # Define quadrature points & weights on [-1,1] interval
    xs = [
        0.0, 0.207784955007898, 0.405845151377397, 0.586087235467691,
        0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813
    ]
    ws = [
        0.209482141084728, 0.204432940075298, 0.190350578064785, 0.169004726639267,
        0.140653259715525, 0.104790010322250, 0.063092092629979, 0.022935322010529
    ]
    xs = np.array(xs[::-1][:-1] + xs)
    ws = np.array(ws[::-1][:-1] + ws)
    xs[:7] *= -1
    
    # Determine the number & size of sub-intervals
    subsize = (bounds[1] - bounds[0]) / nsplit
    subints = [[bounds[0] + subsize * n, bounds[0] + subsize * (n+1)] for n in range(nsplit)]
    
    # Looping over sub-intervals, carry out the integral
    int_out = 0
    for sub in subints:
        
        # Scale quadrature points & weights to desired sub-interval
        sub_xs = (sub[1] - sub[0]) * (xs + 1) / 2 + sub[0]
        sub_ws = (sub[1] - sub[0]) * ws / 2
        
        # Evaluate & weight function at quadrature points along the sub-interval
        int_out += sub_ws @ fun(sub_xs, **fun_args)
    
    return int_out

def EULER(X0, bounds, dt, diffEq, *diffArgs):
    """A simple Euler integrator"""
    
    N = _get_size(bounds, dt)
    isScalar = isinstance(X0,(int,float))
    
    # Intialize state scalar or vector
    X = np.zeros(N) if isScalar else np.zeros((N, len(X0)))
    X[0] = X0

    # Looping over an independent variable, perform Euler integration
    for i in range(N-1):
        t = bounds[0] + i * dt

        # Compute array of slopes
        X[i+1] = X[i] + diffEq(t+dt, X[i], *diffArgs) * dt
        
    return X.tolist()

def RK4(X0, bounds, dt, diffEq, *diffArgs):
    """A standard 4th order Runge-Kutta integrator"""
    
    N = _get_size(bounds, dt)
    isScalar = isinstance(X0,(int,float))

    # Define integrator parameters
    a = np.array([1/2, 1/2, 1])
    b = np.array([1/6, 1/3, 1/3, 1/6])
    c = a

    # Intialize array of state scalars or vectors
    X = np.zeros(N) if isScalar else np.zeros((N, len(X0)))
    X[0] = X0

    # Looping over an independent variable, perform RK4 integration
    for i in range(N-1):
        t = bounds[0]+i*dt

        # Compute array of slopes
        ks = np.zeros(4) if isScalar else np.zeros((4, len(X0)))
        ks[0] = diffEq(t, X[i], *diffArgs)
        for j in range(3):
            if j==0: dXk = sum(ks[:j+1].squeeze() * a[:j+1])
            else: dXk = ks[:j+1].squeeze().T@a[:j+1]
            ks[j+1] = diffEq(t + c[j] * dt, X[i] + dXk * dt, *diffArgs)
        
        # Compute next step using slopes & parameters
        X[i+1] = X[i] + b @ ks * dt

    return X.tolist()

def RK45(X0, bounds, dt, diffEq, *diffArgs):
    """The Dormand-Prince Runge-Kutta integrator"""
    
    N = _get_size(bounds, dt)
    isScalar = isinstance(X0,(int,float))

    # Define integrator parameters
    a = [
        np.array([1/5]),
        np.array([3/40, 9/40]),
        np.array([44/45, -56/15, 32/9]),
        np.array([19372/6561, -25360/2187, 64448/6571, -212/729]),
        np.array([9017/3168, -355/33, 46732/5247, 49/176, -5103/18656]),
        np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84])
    ]
    b = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])
    c = np.array([1/5, 3/10, 4/5, 8/9, 1, 1])

    # Intialize state scalar or vector
    X = np.zeros(N) if isScalar else np.zeros((N,len(X0)))
    X[0] = X0

    # Looping over time, perform RK4 integration
    for i in range(N-1):
        t = bounds[0] + i * dt

        # Compute array of slopes
        ks = np.zeros(7) if isScalar else np.zeros((7, len(X0)))
        ks[0] = diffEq(t, X[i], *diffArgs)
        for j in range(6):
            if j==0: dXk = sum(ks[:j+1].squeeze() * a[j])
            else: dXk = ks[:j+1].squeeze().T @ a[j]
            ks[j+1] = diffEq(t + c[j] * dt, X[i] + dXk * dt, *diffArgs)
        
        # Compute next step using slopes & parameters
        X[i+1] = X[i] + b@ks*dt
        
    return X.tolist()
