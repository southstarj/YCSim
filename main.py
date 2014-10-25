import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
from decimal import *

def calcProp(steam, p):
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
    drhoUw = rhow*steam.diffProp(steam.waterIntEnergy, p)+Uw*drhow
    drhoUg = rhog*steam.diffProp(steam.steamIntEnergy, p)+Ug*drhog
#    print Uw, Ug, drhoUw, drhoUg
    return [rhow, rhog, drhow, drhog,\
            Hw, Hg, Uw, Ug, drhoUw, drhoUg];

def Transmissibility(Sg):
    Tw = 1990; Tg = 10;
    return Tw*(1-Sg) + Tg*Sg;

class Reservoir:
    def __init__(self, n, Vp):
        self._size = n;
        self._poreVol = Vp;

    def Size(self):
        return self._size;
    def checkBound(self, i):
        return (i >= 0) and (i < self._size);
    def poreVolumeList(self):
        return _poreVol;
    def poreVolume(self, i):
        if (self.checkBound(i)):
            return self._poreVol[i];
        else:
            print 'Coordinate out of bound'
            return -1;
n = 2;
Vp = [1000, 1000];
reservoir = Reservoir(n, Vp);
p0 = [600, 500];
Sg0 = [0.0, 1.0];
dt = 1;

steam = prop.steamProp("saturated_steam.org");
getcontext().prec = 5;

"""
dSg=[]
dp=[]
compressibility=[]
steamSaturation=[]
pressure=[]
residual=[]


print '\\begin{table}[H]\n\\centering'
print '\\begin{tabular}{ r r r r r r r r r }'
print '$J$ & $\\delta S_g$ & $\\delta p$ & $\\hat{a}_{22}$'
print '& $S_g$ & $p$ & $\\hat{R}_e$ \\\\'
"""
for J in range(1):
    """
    print '\\centerline{Table',Jvalues.index(J),'}'
    print '\\centerline{Water Injection into Saturated Steam}'
    print '\\vspace{20pt}\n'

    print '\\begin{tabular}{ l l }'
    print '    Injection pressure: &', pi, 'psi \\\\'
    print '    Cell pressure: &', p0, 'psi \\\\'
    print '    Cell gas saturation: &', Sg0, '\\\\'
    print '    Injection enthalpy: &', Hw,\
          'Btu/lb.(Saturated water at 60psi) \\\\'
    print '    Injectivity: &', J, '$lb/ft^3\\cdot psi$ \\\\'
    print '\\end{tabular}'
    print '\\vspace{20pt}\n'

    print '\\begin{table}[H]\n\\centering'
    print '\\begin{tabular}{ r r r r r r r }'
    print 'Itn & $\\delta S_g$ & $\\delta p$ & $\\hat{a}_{22}$'
    print '& $S_g$ & $p$ & $\\hat{R}_e$ \\\\'
    """
    # Prototype for Newton iteration
    Sg = np.array(Sg0);
    p = np.array(p0);
    for iter in range(3):
        A = generateJacobian(reservoir, steam, dt, p, Sg);
        RHS = generateRHS(reservoir, steam, dt, p, Sg);
        x, comp, Rbar = linearSolver(reservoir, steam, dt, A, RHS);
        dp = x[0:1];
        dSg = x[2:3];
        p = p + dp;
        Sg = Sg + dSg;
        print iter, '&',\
              Decimal(dp).normalize(),'&',\
              Decimal(dSg).normalize(),'&',\
              Decimal(comp).normalize(),'&',\
              Decimal(Sg).normalize(),'&',\
              Decimal(p).normalize(),'&',\
              Decimal(Rebar).normalize(),'\\\\'

def generateJacobian(reservoir, steam, dt, p, Sg):
    Vp = reservoir.poreVolumeList();
    rhow = [steam.waterDensity(p[0]), steam.waterDensity(p[1])];
    rhog = [steam.steamDensity(p[0]), steam.steamDensity(p[1])];
    Uw = [steam.waterIntEnergy(p[0]), steam.waterIntEnergy(p[1])];
    Ug = [steam.steamIntEnergy(p[0]), steam.steamIntEnergy(p[1])];
    Hw = [steam.waterEnthalpy(p[0]), steam.waterEnthalpy(p[1])];
    Hg = [steam.steamEnthalpy(p[0]), steam.steamEnthalpy(p[1])];
    T = 1990*(1-Sg[0])+10*Sg[0];
    TH = 1990*(1-Sg[0])*Hw[0]+10*Sg[0]*Hg[0];

    a11 = np.diag([Vp[0]/dt * (rhog[0] - rhow[0]),\
                   Vp[1]/dt * (rhog[1] - rhow[1])]);
    a21 = np.diag([Vp[0]/dt * (rhog[0]*Ug[0] - rhow[0]*Uw[0]),\
                   Vp[1]/dt * (rhog[1]*Ug[1] - rhow[1]*Uw[1])]);


"""
    print '\\end{tabular}\n\\end{table}'
    print '\\newpage'
"""
