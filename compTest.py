import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
import Reservoir
import Fluid
from Simulator import *

n = 1;
nv = 2;
debugstep = -1;
# geological settings
poreVol = [1000.0 for i in range(n)];
perm = [0.5 for i in range(n)];
deltax = 100.0;
Area = 10.0;
reservoir = Reservoir.Reservoir(n, poreVol, perm, deltax, Area);
# initial distribution
p0 = [500 for i in range(n)];  p0[0] = 300;
Sg0 = [1.0 for i in range(n)]; Sg0[0] = 1.0;
Wellnum = 0; qT = -5.0; pWater = .2563; pInj = 300;            # water injection

rhoB = reservoir.getFluid().waterDensity(pInj);
HB = reservoir.getFluid().waterEnthalpy(pWater);
muB = reservoir.getFluid().waterViscosity(pWater);
T = reservoir.Permeability(Wellnum)*reservoir.GetSectionA()\
    /(reservoir.Deltax(Wellnum)*muB);


x0 = np.array(Sg0 + p0);

dt = 1;                           # time step

np.set_printoptions(precision=5);

for p in [190, 195, 198, 199]:
    x = np.array([1.0, p]);
    x0 = x;
    RHS = GenerateRHS(reservoir, dt, x, x0);
    #RHS = BoundaryCond_Rate(reservoir, RHS, qT, pWater, pInj, Wellnum);
    A = GenerateJacobian(reservoir, dt, x);
    A, RHS = BoundaryCond_Pres(reservoir, A, RHS, pWater, pInj,\
             x[reservoir.Size()+Wellnum], Wellnum);
    dx, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
    print 'p =', p
    print '  Injectivity =', T*(pInj - p)/poreVol[0]
    print '  A ='
    print A
    print '  RHS =', RHS
    print '  dx =', dx[0], dx[1]
    print '  comp =', comp
    print '  Rebar =', Rebar