import numpy as np
import scipy as sp
import matplotlib as plt

pi = 600;      # inject pressure: psia
p0 = 500;      # cell pressure: psia
Sg0 = 1;       # init saturation
Hw = 262.09;   # inject water enthalpy: btu/lb
J = 0.005;     # normalized injectivity

crho = 0.062428*0.178107607;    # kg/m3 to lb/RB
rhow = 810.96*crho; # water density: kg/m3
rhog = 17.259*crho; # steam density: kg/m3
drhow = (810.80-811.13)/2*crho; # drhow/dp
drhog = (17.294-17.224)/2*crho; # drhog/dp

cU = 0.429923;      # kJ/kg to Btu/lb
Uw = 1041.37*cU;    # water internal energy: kJ/kg
Ug = 2603.12*cU;    # steam internal energy: kJ/kg
drhoUw = (810.80*1041.91-811.13*1040.83)/2*crho*cU;  # d(rhowUw)/dp
drhoUg = (17.294*2603.11-17.224*2603.13)/2*crho*cU;  # d(rhogUg)/dp
print drhoUw, drhoUg

Sg = Sg0;
p = p0;

for iter in range(5):
    Rw = (rhow*(1-Sg) + rhog*Sg) - J*(pi-p);
    Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - J*Hw*(pi-p);
    RHS = -np.array([Rw, Re]).transpose();
    a11 = rhog-rhow;
    a21 = rhog*Ug-rhow*Uw;
    a12 = (1-Sg)*drhow + Sg*drhog + J;
    a22 = (1-Sg)*drhoUw + Sg*drhoUg + J*Hw;
    comp = a22 - a21*a12/a11;
    Rebar = Re - a21/a11*Rw;
#    print a11, a21, a12, a22
    A = np.array([[a11, a12], [a21, a22]]);
    print A
    x = np.linalg.solve(A, RHS);
    Sg = Sg0 + x[0];
    p = p0 + x[1];
    print x[0], x[1], comp, Sg, p, Rebar;
