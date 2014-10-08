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
    print 'rhow =', rhow
    print 'rhog =', rhog
    print 'drhow/dp =', drhow
    print 'drhog/dp =', drhog
    cU = 0.429923;      # kJ/kg to Btu/lb
    Hw = steam.waterEnthalpy(p)
    Hg = steam.steamEnthalpy(p)
    Uw = steam.waterIntEnergy(p)
    Ug = steam.steamIntEnergy(p)
    drhoUw = rhow*steam.diffProp(steam.waterIntEnergy, p)+Uw*drhow
    drhoUg = rhog*steam.diffProp(steam.steamIntEnergy, p)+Ug*drhog
    print 'Hw =', Hw
    print 'Hg =', Hg
    print 'Uw =', Uw
    print 'Ug =', Ug
    print 'drhowUw/dp =', drhoUw
    print 'drhogUg/dp =', drhoUg
#    return [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg];

steam = prop.steamProp("saturated_steam.org");
p = 500;
print 'Steam properties at', p, 'psi'
calcProp(steam, p)
