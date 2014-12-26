import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
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
    #print 'in QTerm:', 'p =', p, 'Sg =', Sg
    steam = reservoir.getFluid();
    rhow = steam.waterDensity(p);
    rhog = steam.steamDensity(p);
    Uw = steam.waterIntEnergy(p);
    Ug = steam.steamIntEnergy(p);
    # calculate function value
    if diffVar == 0:
        Qw = []; Qe = [];
        for i in range(n):
            Qw.append(Vp[i]/dt*(rhow[i]*(1-Sg[i])+rhog[i]*Sg[i]));
            Qe.append(Vp[i]/dt*(rhow[i]*(1-Sg[i])*Uw[i]+rhog[i]*Sg[i]*Ug[i]));
        return (Qw, Qe);
    # calculate differential
    dQw = []; dQe = [];
    if diffVar == 2:        # differential p
        drhow = steam.diffProp(steam.waterDensity, p);
        drhog = steam.diffProp(steam.steamDensity, p);
        dUw = steam.diffProp(steam.waterIntEnergy, p);
        dUg = steam.diffProp(steam.steamIntEnergy, p);
        drhoUw = [(x*dy+y*dx) for x,y,dx,dy in zip(rhow,Uw,drhow,dUw)];
        drhoUg = [(x*dy+y*dx) for x,y,dx,dy in zip(rhog,Ug,drhog,dUg)];
        for i in range(n):
            dQw.append(Vp[i]/dt*(drhow[i]*(1-Sg[i])+drhog[i]*Sg[i]));
            dQe.append(Vp[i]/dt*(drhoUw[i]*(1-Sg[i])+drhoUg[i]*Sg[i]));
    if diffVar == 1:        # differential Sg
        for i in range(n):
            dQw.append(Vp[i]/dt*(rhog[i] - rhow[i]));
            dQe.append(Vp[i]/dt*(rhog[i]*Ug[i] - rhow[i]*Uw[i]));
    return (dQw, dQe);

def AccumulationTerm(reservoir, dt, p, Sg, p0, Sg0):
    Qw1, Qe1 = QTerm(reservoir, dt, p, Sg);
    Qw0, Qe0 = QTerm(reservoir, dt, p0, Sg0);
    #print 'Qwn, Qw0, Qen, Qe0 =', Qw1, Qw0, Qe1, Qe0

    return (np.array(Qw1) - np.array(Qw0), np.array(Qe1) - np.array(Qe0));

def GenerateRHS(reservoir, dt, x, x0):
    n = reservoir.Size();
    Sg = x[0:n];
    p = x[n:2*n];
    Sg0 = x0[0:n];
    p0 = x0[n:2*n];

    DQw, DQe = AccumulationTerm(reservoir, dt, p, Sg, p0, Sg0);
    Rw = DQw; Re = DQe;
    #print 'DQw, DQe =', Rw, Re
    #print 'p, p0, Sg, Sg0 =', p, p0, Sg, Sg0
    for i in range(n):
        _conn = reservoir.GetConnection(i);
        for k in _conn:
            T, TH = reservoir.Transmissibility(i, k, x);
            Dp = p[k] - p[i];
            #print 'i, k = (', i, k, '), T =', T, 'Dp =', Dp
            Rw[i] += -np.sum(T)*Dp;
            Re[i] += -np.sum(TH)*Dp;

    RHS = -np.hstack((Rw, Re));
    #print Rw, Re
    return RHS;

def GenerateJacobian(reservoir, dt, x):
    n = reservoir.Size();
    Vp = reservoir.PoreVolumeList();
    Sg = x[0:n];
    p = x[n:2*n];

    dQw_dS, dQe_dS = QTerm(reservoir, dt, p, Sg, 1);
    dQw_dp, dQe_dp = QTerm(reservoir, dt, p, Sg, 2);
    a11 = np.diag(dQw_dS);
    a21 = np.diag(dQe_dS);
    a12 = np.diag(dQw_dp);
    a22 = np.diag(dQe_dp);
    """
    print 'a11 =', a11
    print 'a12 =', a12
    print 'a21 =', a21
    print 'a22 =', a22
    """
    for i in range(n):
        _conn = reservoir.GetConnection(i);
        for k in _conn:
            Dp = p[k] - p[i];
            T, TH = reservoir.Transmissibility(i, k, x);
            dT_dSk, dTH_dSk = reservoir.Transmissibility(k, i, x, 1);
            dT_dpk, dTH_dpk = reservoir.Transmissibility(k, i, x, 2);
            dT_dSi, dTH_dSi = reservoir.Transmissibility(i, k, x, 1);
            dT_dpi, dTH_dpi = reservoir.Transmissibility(i, k, x, 2);
            """
            print 'i, k = (', i, k, '):', 'Dp =', Dp
            print '  dT_dSk', dT_dSk, 'dTH_dSk', dTH_dSk
            print '  dT_dpk', dT_dpk, 'dTH_dpk', dTH_dpk
            print '  dT_dSi', dT_dSi, 'dTH_dSi', dTH_dSi
            print '  dT_dpi', dT_dpi, 'dTH_dpi', dTH_dpi
            """
            a11[i, k] = -np.sum(dT_dSk)*Dp;
            a11[i, i] -= np.sum(dT_dSi)*Dp;
            a21[i, k] = -np.sum(dTH_dSk)*Dp;
            a21[i, i] -= np.sum(dTH_dSi)*Dp;
            a12[i, k] = -np.sum(T) - np.sum(dT_dpk)*Dp;
            a12[i, i] += np.sum(T) - np.sum(dT_dpi)*Dp;
            a22[i, k] = -np.sum(TH) - np.sum(dTH_dpk)*Dp;
            a22[i, i] += np.sum(TH) - np.sum(dTH_dpi)*Dp;

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
    """
    print 'a11 ='
    print a11
    print 'a12 ='
    print a12
    print 'a21 ='
    print a21
    print 'a22 ='
    print a22
    """
    print 'factors ='
    print np.dot(_factor, a12), np.dot(_factor, Rw)
    print 'comp ='
    print comp
    print 'Rebar ='
    print Rebar
    return (dx, comp, Rebar);

def BoundaryCond_Rate(reservoir, RHS, qT, pB):
    rhoB = reservoir.getFluid().waterDensity(pB);
    HB = reservoir.getFluid().waterEnthalpy(pB);
    #print 'rhoB, HB, qT =', rhoB, HB, qT
    #print RHS
    RHS[0] += -qT*rhoB;
    RHS[reservoir.Size()] += -qT*rhoB*HB;
    return RHS;
