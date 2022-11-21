### integrator_tests.py
### Morgan Henderson, July 2022
### Tests various integrators used in the hysteresis.ipynb notebook

from numpy import array, sin, log, exp, pi, einsum
from integrators import EULER, RK4, RK45, KRON15
e = exp(1)

# Define lambda functions to alert user of test status
warn = lambda s: print('WARNING: %s integration not functioning properly.'%s)
success = lambda s: print('SUCCESS: %s integration functioning properly.'%s)

# Define test functions for quadrature integration
quadFuns = [lambda x: 3*x**2, lambda x: sin(x), lambda x: log(x)/x]
quadBounds = [(0,1), (0,pi), (1,e)]
quadAns = [1, 2, 0.5]

# Define scalar & vector differential equations to test integrators
def intTest1(t,x):
    return -exp(-t)
def intTest2(t,x,w=1):
    return array([x[1],-x[0]*w**2])
tInt = array(range(10001))/1000
x0Int1 = 1
x0Int2 = array([1,0])
intAns1 = exp(-tInt).T
intAns2 = array([sin(tInt+pi/2), -sin(tInt)]).T

# Test the KRON15 quadrature integrator on various functions
kronAcc = 3
kronAns = [KRON15(quadFuns[t],quadBounds[t]) for t in range(3)]
kronTests = [round(kronAns[t],kronAcc)==quadAns[t] for t in range(3)]
if not all(kronTests): warn('KRON15 quadrature')
else: success('KRON15 quadrature')

# Test the EULER integrator against a critical error value
eulerETol = 0.05
eulerAns = [EULER(x0Int1,tInt[[0,-1]],tInt[1]-tInt[0],intTest1),\
              EULER(x0Int2,tInt[[0,-1]],tInt[1]-tInt[0],intTest2)]
eulerErr = [eulerAns[0]-intAns1, eulerAns[1]-intAns2]
eulerTests = [(eulerErr[0]@eulerErr[0]/len(tInt))**.5,\
             (einsum('ij,ij->j',eulerErr[1],eulerErr[1])/len(tInt))**.5]
if not eulerTests[0]<eulerETol: warn('EULER scalar')
else: success('EULER scalar')
if not all(eulerTests[1]<eulerETol): warn('EULER vector')
else: success('EULER vector')

# Test the RK4 & RK45 integrators with scalar & vector states
rkETol = [1e-3,1e-2]
rk4Ans = [RK4(x0Int1,tInt[[0,-1]],tInt[1]-tInt[0],intTest1),\
              RK4(x0Int2,tInt[[0,-1]],tInt[1]-tInt[0],intTest2)]
rk4Err = [rk4Ans[0]-intAns1, rk4Ans[1]-intAns2]
rk4Tests = [(rk4Err[0]@rk4Err[0]/len(tInt))**.5,\
             (einsum('ij,ij->j',rk4Err[1],rk4Err[1])/len(tInt))**.5]
if not rk4Tests[0]<rkETol[0]: warn('RK4 scalar')
else: success('RK4 scalar')
if not all(rk4Tests[1]<rkETol[1]): warn('RK4 vector')
else: success('RK4 vector')
rk45Ans = [RK45(x0Int1,tInt[[0,-1]],tInt[1]-tInt[0],intTest1),\
              RK45(x0Int2,tInt[[0,-1]],tInt[1]-tInt[0],intTest2)]
rk45Err = [rk45Ans[0]-intAns1, rk45Ans[1]-intAns2]
rk45Tests = [(rk45Err[0]@rk45Err[0]/len(tInt))**.5,\
             (einsum('ij,ij->j',rk45Err[1],rk45Err[1])/len(tInt))**.5]
if not rk45Tests[0]<rkETol[0]: warn('RK45 scalar')
else: success('RK45 scalar')
if not all(rk45Tests[1]<rkETol[1]): warn('RK45 vector')
else: success('RK45 vector')