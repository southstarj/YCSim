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
p0 = [500 for i in range(n)];  p0[0] = 200;
Sg0 = [1.0 for i in range(n)]; Sg0[0] = 1.0;
Wellnum = 0; qT = -5.0; pWater = .2563; pInj = 300;            # water injection
rhoB = reservoir.getFluid().waterDensity(pInj);
HB = reservoir.getFluid().waterEnthalpy(pWater);
muB = reservoir.getFluid().waterViscosity(pWater);
J = reservoir.Permeability(Wellnum)*reservoir.GetSectionA()\
    /(reservoir.Deltax(Wellnum)*muB);
print 'Injectivity =', J*rhoB, J*rhoB*HB
"""
print 'Water Injection:', qT, 'ft^3/h'
print rhoB, 'lbm/ft^3 (', pWater, 'psi)'
print HB, 'Btu/lbm (', pInj, 'psi)\n'
"""
x0 = np.array(Sg0 + p0);
x = x0;
dt = 4;                           # time step

np.set_printoptions(precision=5);

#Legend = [];
"""
fig = plt.figure(figsize=(16, 9), dpi = 80);
AxSaturation = fig.add_subplot(412)
AxGridblock = fig.add_subplot(411, sharex=AxSaturation)
AxPressure = fig.add_subplot(413, sharex=AxSaturation)
AxTemperature = fig.add_subplot(414, sharex=AxSaturation)
"""
tt = [];
pt = [];
Sgt = [];
for timestep in range(101):
    print 'time step:', timestep
    # Prototype for Newton iteration
    for iter in range(100):
        #print 'iter =', iter
        #T, TH = reservoir.Transmissibility(0, 1, x)

        RHS = GenerateRHS(reservoir, dt, x, x0);
        #RHS = BoundaryCond_Rate(reservoir, RHS, qT, pWater, pInj, Wellnum);
        A = GenerateJacobian(reservoir, dt, x);
        A, RHS = BoundaryCond_Pres(reservoir, A, RHS, pWater, pInj,\
                 x[reservoir.Size()+Wellnum], Wellnum);
        if np.linalg.norm(RHS) < 1e-3:
            print '  iteration number:', iter
            break;
        #print A
        #print RHS
        dx = np.linalg.solve(A, RHS);
        #x0 = x;
        #x = x + dx;
        #print 'iter =', iter, x

        #dx, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
        x = x + dx;
        if timestep == debugstep:
            print 'iter =', iter
            #J = -qT*rhoB*dt/poreVol[0]/(pInj-x[n]); JH = J*HB;
            #J = np.sum(T)*dt/poreVol[0]; JH = J*HB;
            #print 'injectivity =', J, JH
            #print 'Transmissibility =', T, TH
            #print 'Dp =', x[2] - x[3]
            print 'Jacobian ='
            print A #np.array([[A[0][0], A[0][2]],[A[2][0], A[2][2]]])
            print 'RHS =', RHS #np.array([RHS[0], RHS[2]])
            dx, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
            print 'comp =', comp
            print 'Rebar =', Rebar
            print 'dx =', dx #np.array([dx[0], dx[2]])
            #print 'new dx =', dx
            print 'x =', x #np.array([x[0], x[2]])

    if timestep == debugstep:
        print x0
        #print dx
        print x
        break
    Sg = x[0:n];
    p = x[n:(2*n)];
    Tb = reservoir.getFluid().boilingPoint(p);
    print 'Sg =', Sg
    print 'p =', p
    tt.append(timestep);
    pt.append(p[0]);
    Sgt.append(Sg[0]);
    x0 = x;
    if iter == 0:
        print 'Steady State'
        break;
    if iter == 99:
        print 'Not Converge'
        break;
    """
    if timestep%20 == 0:
        AxSaturation.plot(Sg, ':o')
        AxPressure.plot(p, ':o')
        AxTemperature.plot(Tb, ':o')
        Legend.append('t = ' + str(timestep+1))
    """
"""
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
"""
"""
AxGridblock.annotate('$q_T,\;p_{inj},\;H_{w,inj}$',
            xy=(0, 0.5), xycoords='data',
            xytext=(-1, 0.5), textcoords='data',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->")
           )

AxGridblock.annotate('$p_{prod}$',
            xy=(9, 0.5), xycoords='data',
            xytext=(8.5, 0.5), textcoords='data',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->")
           )

AxGridblock.annotate('$A,\Delta x$', xy=(3.5, 0.5), xycoords='data',
            bbox=dict(boxstyle="round", fc="0.8"),
           )

AxGridblock.annotate('$k, V_p$', xy=(4.5, 0.5), xycoords='data',
            bbox=dict(boxstyle="round", fc="0.8"),
           )
"""
#AxGridblock.set_xbound(0, 9)

fig = plt.figure();
AxPressure = fig.add_subplot(211);
#plt.title('$q_T='+str(-qT)+'ft^3/d\;H_{w,inj}=262.25btu/lbm(60psi)$')
plt.title('Injectivity = '+str(J*rhoB/1000*dt)+' '+str(J*rhoB*HB/1000*dt));
AxSaturation = fig.add_subplot(212);
AxPressure.plot(pt)
AxSaturation.plot(Sgt)
AxSaturation.grid(True)
AxSaturation.set_ylabel('Sg')
AxPressure.grid(True)
AxPressure.set_ylabel('p(psi)')
#AxPressure.set_ybound(480, 500)

plt.show()