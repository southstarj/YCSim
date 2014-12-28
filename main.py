import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
import Reservoir
import Fluid
from Simulator import *

n = 10;
nv = 2;
# geological settings
poreVol = [100.0 for i in range(n)];
perm = [2.0 for i in range(n)];
deltax = 10.0;
Area = 10.0;
reservoir = Reservoir.Reservoir(n, poreVol, perm, deltax, Area);
# initial distribution
p0 = [500 for i in range(n)];  p0[0] = 500;
Sg0 = [1.0 for i in range(n)]; Sg0[0] = 1.0;
qT = -1.0; pWater = 600; pInj = 600;                # water injection
rhoB = reservoir.getFluid().waterDensity(pInj);
HB = reservoir.getFluid().waterEnthalpy(pWater);
x0 = np.array(Sg0 + p0);
x = x0;
dt = 10.0;                           # time step

np.set_printoptions(precision=5);

Legend = [];
fig = plt.figure(figsize=(16, 9), dpi = 80);
AxSaturation = fig.add_subplot(412)
AxGridblock = fig.add_subplot(411, sharex=AxSaturation)
AxPressure = fig.add_subplot(413, sharex=AxSaturation)
AxTemperature = fig.add_subplot(414, sharex=AxSaturation)

for timestep in range(201):
    print 'time step:', timestep
    # Prototype for Newton iteration
    for iter in range(100):
        #print 'iter =', iter
        T, TH = reservoir.Transmissibility(0, 1, x)

        RHS = GenerateRHS(reservoir, dt, x, x0);
        RHS = BoundaryCond_Rate(reservoir, RHS, qT, pWater, pInj);
        if np.linalg.norm(RHS) < 1e-3:
            print '  iteration number:', iter
            break;
        A = GenerateJacobian(reservoir, dt, x);
        #print A
        #print RHS
        dx = np.linalg.solve(A, RHS);
        #x0 = x;
        #x = x + dx;
        #print 'iter =', iter, x

        #dx, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
        x = x + dx;
        if timestep == -1:
            print 'iter =', iter
            J = qT*rhoB*dt/poreVol[0]/(pInj-x[2]); JH = J*HB;
            print 'injectivity =', J, JH
            print 'Transmissibility =', T, TH
            print 'Dp =', x[2] - x[3]
            print 'Jacobian ='
            print np.array([[A[0][0], A[0][2]],[A[2][0], A[2][2]]])
            print 'RHS =', np.array([RHS[0], RHS[2]])
            dx, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
            print 'comp =', comp[0][0]
            print 'Rebar =', Rebar[0]
            print 'dx =', np.array([dx[0], dx[2]])
            #print 'new dx =', dx
            print 'x =', np.array([x[0], x[2]])

    if timestep == -1:
        print x0
        #print dx
        print x
        break
    Sg = x[0:n];
    p = x[n:(2*n)];
    Tb = reservoir.getFluid().boilingPoint(p);
    print 'Sg =', Sg
    print 'p =', p
    x0 = x;
    if iter == 0:
        print 'Steady State'
        break;
    if iter == 99:
        print 'Not Converge'
        break;
    
    if timestep%20 == 0:
        AxSaturation.plot(Sg)
        AxPressure.plot(p)
        AxTemperature.plot(Tb)
        Legend.append('t = ' + str(timestep+1))

AxSaturation.grid(True)
AxSaturation.set_ylabel('Sg')
AxPressure.grid(True)
AxPressure.set_ylabel('p(psi)')
AxTemperature.grid(True)
AxTemperature.set_ylabel('Tb(F)')
AxTemperature.set_xlabel('Grid block')
AxTemperature.legend(Legend, fontsize = 8)
AxGridblock.set_yticks([])
AxGridblock.grid(True)
AxGridblock.set_xbound(0, 9)

#plt.show()