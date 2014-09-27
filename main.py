import numpy as np
import scipy as sp
import matplotlib as plt
import steamProp as prop

def calcProp(steam, p):
    crho = 0.062428*0.178107607;    # kg/m3 to lb/RB
    rhow = steam.waterDensity(p)
    rhog = steam.steamDensity(p)
    drhow = steam.diffProp(steam.waterDensity, p)
    drhog = steam.diffProp(steam.steamDensity, p)
#    print rhow, rhog, drhow, drhog
    cU = 0.429923;      # kJ/kg to Btu/lb
    Uw = steam.waterIntEnergy(p)
    Ug = steam.steamIntEnergy(p)
    drhoUw = rhow*steam.diffProp(steam.waterIntEnergy, p)+Uw*drhow
    drhoUg = rhog*steam.diffProp(steam.steamIntEnergy, p)+Ug*drhog
#    print Uw, Ug, drhoUw, drhoUg
    return [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg];

pi = 600;      # inject pressure: psia
p0 = 500;      # cell pressure: psia
Sg0 = 1.0;     # init saturation
Hw = 471.6;   # inject water enthalpy: btu/lb
J = 10.0;     # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");

Sg = Sg0;
p = p0;
[rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
Rw0 = rhow*(1-Sg) + rhog*Sg;
Re0 = (rhow*(1-Sg)*Uw + rhog*Sg*Ug);
v = calcProp(steam, 595.03);
print 'Uw(595.03psi) =',v[4]

for iter in range(30):
#    print 'Iteration #', iter
    Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
    Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
    RHS = -np.array([Rw, Re]).transpose();
#    print 'RHS = ', RHS
    a11 = rhog-rhow;
    a21 = rhog*Ug-rhow*Uw;
    a12 = (1-Sg)*drhow + Sg*drhog + J;
    a22 = (1-Sg)*drhoUw + Sg*drhoUg + J*Hw;
    comp = a22 - a21/a11*a12;
    Rebar = Re - a21/a11*Rw;
#    print a11, a21, a12, a22
    A = np.array([[a11, a12], [a21, a22]]);
#    print A
    x = np.linalg.solve(A, RHS);
    Sg = Sg + x[0];
    p = p + x[1];
    """
    print '    dSg =', x[0]
    print '    dp =', x[1]
    print '    Compressibility =', comp
    print '  Sg =', Sg
    print '  p =', p
    print '    modified Re =', Rebar
    """
    print iter, '\t', x[0],'\t', x[1],'\t', comp,'\t', Sg,'\t', p,'\t', Rebar
    [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
