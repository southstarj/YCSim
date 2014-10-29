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
    dUw = steam.diffProp(steam.waterIntEnergy, p);
    dUg = steam.diffProp(steam.steamIntEnergy, p);
    drhoUw = [(x*dy+y*dx) for x,y,dx,dy in zip(rhow,Uw,drhow,dUw)];
    drhoUg = [(x*dy+y*dx) for x,y,dx,dy in zip(rhog,Ug,drhog,dUg)];
#    print Uw, Ug, drhoUw, drhoUg
    return (rhow, rhog, drhow, drhog,\
            Hw, Hg, Uw, Ug, drhoUw, drhoUg);

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
        return self._poreVol;
    def poreVolume(self, i):
        if (self.checkBound(i)):
            return self._poreVol[i];
        else:
            print 'Coordinate out of bound'
            return -1;


def accumulationTerm(reservoir, steam, dt, p, Sg, p0, Sg0):
    Vp = reservoir.poreVolumeList();
    (rhow, rhog, drhow, drhog, Hw, Hg, Uw, Ug, drhoUw, drhoUg) =\
        calcProp(steam, p);
    (rhow0, rhog0, drhow, drhog, Hw0, Hg0, Uw0, Ug0, drhoUw, drhoUg) =\
        calcProp(steam, p0);
    Qw = []; Qe = [];

    for i in range(reservoir.Size()):
        Qw.append(Vp[i]/dt*((rhow[i]*(1-Sg[i])+rhog[i]*Sg[i]) -\
                            (rhow0[i]*(1-Sg0[i])+rhog0[i]*Sg0[i])));
        Qe.append(Vp[i]/dt*((rhow[i]*(1-Sg[i])*Uw[i]+rhog[i]*Sg[i]*Ug[i])-\
                     (rhow0[i]*(1-Sg0[i])*Uw0[i]+rhog0[i]*Sg0[i]*Ug0[i])));
    return (Qw, Qe);

def generateRHS(reservoir, steam, dt, p, Sg, p0, Sg0):
    (rhow, rhog, drhow, drhog, Hw, Hg, Uw, Ug, drhoUw, drhoUg) =\
        calcProp(steam, p);
    T = 1990*(1-Sg[0])+10*Sg[0];
    TH = 1990*(1-Sg[0])*Hw[0]+10*Sg[0]*Hg[0];
    Qw, Qe = accumulationTerm(reservoir, steam, dt, p, Sg, p0, Sg0);

    Dp = p[0]-p[1];
    Rw0 = T*Dp + Qw[0];
    Rw1 = -T*Dp+ Qw[1];
    Re0 = TH*Dp+ Qe[0];
    Re1 =-TH*Dp+ Qe[1];
    RHS = -np.array([Rw0, Rw1, Re0, Re1]);
    return RHS

def generateJacobian(reservoir, steam, dt, p, Sg):
    Vp = reservoir.poreVolumeList();
    (rhow, rhog, drhow, drhog, Hw, Hg, Uw, Ug, drhoUw, drhoUg) =\
        calcProp(steam, p);
    T = 1990*(1-Sg[0])+10*Sg[0];
    TH = 1990*(1-Sg[0])*Hw[0]+10*Sg[0]*Hg[0];

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
    p0 = np.array(p0);
    Sg0 = np.array(Sg0);
    Sg = Sg0;
    p = p0;
    for iter in range(1):
        A = generateJacobian(reservoir, steam, dt, p, Sg);
        RHS = generateRHS(reservoir, steam, dt, p, Sg, p0, Sg0);
        print A
        print RHS
        x = np.linalg.solve(A, RHS);
        print x

        #x, comp, Rbar = linearSolver(reservoir, steam, dt, A, RHS);
        dSg = x[0:2];
        dp = x[2:4];
        Sg = Sg + dSg;
        p = p + dp;
        """
        print iter, '&',\
              Decimal(dp).normalize(),'&',\
              Decimal(dSg).normalize(),'&',\
              Decimal(comp).normalize(),'&',\
              Decimal(Sg).normalize(),'&',\
              Decimal(p).normalize(),'&',\
              Decimal(Rebar).normalize(),'\\\\'
        """

"""
    print '\\end{tabular}\n\\end{table}'
    print '\\newpage'
"""
