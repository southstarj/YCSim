import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
from decimal import *
import Reservoir
import Fluid

def CalcProp(steam, p):
    crho = 0.062428*0.178107607;    # kg/m3 to lb/RB
    rhow = steam.waterDensity(p)
    rhog = steam.steamDensity(p)
    drhow = steam.diffProp(steam.waterDensity, p)
    drhog = steam.diffProp(steam.steamDensity, p)
#    print rhow, rhog, drhow, drhog
    cU = 0.429923;      # kJ/kg to Btu/lb
    Hw = steam.waterEnthalpy(p)
    Hg = steam.steamEnthalpy(p)
    Uw = steam.waterIntEnergy(p)
    Ug = steam.steamIntEnergy(p)
    dUw = steam.diffProp(steam.waterIntEnergy, p);
    dUg = steam.diffProp(steam.steamIntEnergy, p);
    drhoUw = [(x*dy+y*dx) for x,y,dx,dy in zip(rhow,Uw,drhow,dUw)];
    drhoUg = [(x*dy+y*dx) for x,y,dx,dy in zip(rhog,Ug,drhog,dUg)];
#    print Uw, Ug, drhoUw, drhoUg
    return (rhow, rhog, drhow, drhog,\
            Hw, Hg, Uw, Ug, drhoUw, drhoUg);

def QTerm(reservoir, dt, p, Sg, diffVar=0):
    n = reservoir.Size();
    Vp = reservoir.PoreVolumeList();
    steam = reservoir.getFluid();
    rhow = steam.waterDensity(p);
    rhog = steam.steamDensity(p);
    Uw = steam.waterIntEnergy(p);
    Ug = steam.steamIntEnergy(p);
    # calculate function value
    if diffVar == 0:
        Qw = []; Qe = [];
        for i in range():
            Qw.append(Vp[i]/dt*(rhow[i]*(1-Sg[i])+rhog[i]*Sg[i]));
            Qe.append(Vp[i]/dt*(rhow[i]*(1-Sg[i])*Uw[i]+rhog[i]*Sg[i]*Ug[i]));
        return (Qw, Qe);
    # calculate differential
    dQw = []; dQe = [];
    if diffVar == 1:
        drhow = steam.diffProp(steam.waterDensity, p);
        drhog = steam.diffProp(steam.steamDensity, p);
        drhoUw = [(x*dy+y*dx) for x,y,dx,dy in zip(rhow,Uw,drhow,dUw)];
        drhoUg = [(x*dy+y*dx) for x,y,dx,dy in zip(rhog,Ug,drhog,dUg)];
        for i in range():
            dQw.append(Vp[i]/dt*(drhow[i]*(1-Sg[i])+drhog[i]*Sg[i]));
            dQe.append(Vp[i]/dt*(drhoUw[i]*(1-Sg[i])+drhoUg[i]*Sg[i]));
    if diffVar == 2:
        for i in range():
            dQw.append(Vp[i]/dt*(rhog[i] - rhow[i]));
            dQe.append(Vp[i]/dt*(rhog[i]*Ug[i] - rhow[i]*Uw[i]));
    return (dQw, dQe);

def AccumulationTerm(reservoir, dt, p, Sg, p0, Sg0):
    Qw1, Qe1 = QTerm(reservoir, p, Sg);
    Qw0, Qe0 = QTerm(reservoir, p0, Sg0);

    return (np.array(Qw1) - np.array(Qw0), np.array(Qe1) - np.array(Qe0));

def GenerateRHS(reservoir, dt, x, x0):
    n = reservoir.Size();
    Sg = x[0:n];
    p = x[n:2*n];
    Sg0 = x0[0:n];
    p0 = x0[n:2*n];

    Qw, Qe = AccumulationTerm(reservoir, dt, p, Sg, p0, Sg);
    Rw = Qw; Re = Qe;
    for i in range(n):
        _conn = reservoir.GetConnection(i);
        for k in _conn:
            T, TH = reservoir.Transmissibility(i, k, x);
            Dp = p[k] - p[i];
            Rw[i] += -np.sum(T)*Dp;
            Re[i] += -np.sum(TH)*Dp;

    RHS = -np.hstack((Rw, Re));
    return RHS;

def GenerateJacobian(reservoir, dt, x):
    n = reservoir.Size();
    Vp = reservoir.PoreVolumeList();
    Sg = x[0:n];
    p = x[n:2*n];

    a11 = np.zeros((n, n));
    a12 = np.zeros((n, n));
    a21 = np.zeros((n, n));
    a22 = np.zeros((n, n));
    for i in range(n):
        _conn = reservoir.GetConnection(i);
        for k in _conn:
        # working here

    a11 = np.diag([Vp[0]/dt * (rhog[0] - rhow[0]),\
                   Vp[1]/dt * (rhog[1] - rhow[1])]);
    a21 = np.diag([Vp[0]/dt * (rhog[0]*Ug[0] - rhow[0]*Uw[0]),\
                   Vp[1]/dt * (rhog[1]*Ug[1] - rhow[1]*Uw[1])]);
    a12 = np.array([[T + Vp[0]/dt*((1-Sg[0])*drhow[0]\
                                     +Sg[0] *drhog[0]), -T],\
                    [-T, T + Vp[1]/dt*((1-Sg[1])*drhow[1]\
                                     +Sg[1] *drhog[1])]]);
    a22 = np.array([[TH + Vp[0]/dt*((1-Sg[0])*drhoUw[0]\
                                      +Sg[0] *drhoUg[0]), -TH],\
                    [-TH, TH + Vp[1]/dt*((1-Sg[1])*drhoUw[1]\
                                           +Sg[1] *drhoUg[1])]]);
    A = np.vstack([np.hstack([a11,a12]),\
                   np.hstack([a21,a22])]);
    return A;

def LinearSolver(reservoir, dt, A, RHS):
    a11 = A[0:2, 0:2]; a12 = A[0:2, 2:4];
    a21 = A[2:4, 0:2]; a22 = A[2:4, 2:4];
    Rw = -RHS[0:2]; Re = -RHS[2:4];
    _factor = np.dot(a21, np.linalg.inv(a11));
    comp = a22 - np.dot(_factor, a12);
    Rebar = Re - np.dot(_factor, Rw);
    dp = np.linalg.solve(comp, -Rebar);
    dSg = np.linalg.solve(a11, -Rw - np.dot(a12, dp));
    dx = np.concatenate((dSg, dp));
    print 'a11 ='
    print a11
    print 'a12 ='
    print a12
    print 'a21 ='
    print a21
    print 'a22 ='
    print a22
    print 'factors ='
    print np.dot(_factor, a12), np.dot(_factor, Rw)
    print 'comp ='
    print comp
    print 'Rebar ='
    print Rebar
    return (dx, comp, Rebar);

n = 10;
nv = 2;
Vp = [1000 for i in range(10)];
k = [209 for i in range(10)];
reservoir = Reservoir.Reservoir();
p0 = [600, 500, 500, 500, 500, 500, 500, 500, 500, 400];
Sg0 = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
x0 = np.array(Sg0 + p0);
x = x0;
dt = 1;

np.set_printoptions(precision=5);

for timestep in range(1):
    # Prototype for Newton iteration
    for iter in range(1):
        A = GenerateJacobian(reservoir, dt, x);
        RHS = GenerateRHS(reservoir, dt, x, x0);
        #print A
        #print RHS
        #x = np.linalg.solve(A, RHS);
        #print x

        x, comp, Rebar = LinearSolver(reservoir, dt, A, RHS);
        dSg = x[0:n];
        dp = x[n:(2*n)];
        Sg = Sg + dSg;
        p = p + dp;
