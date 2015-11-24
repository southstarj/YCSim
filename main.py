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
    Uw = steam.waterIntEnergy(p)
    Ug = steam.steamIntEnergy(p)
    drhoUw = rhow*steam.diffProp(steam.waterIntEnergy, p)+Uw*drhow
    drhoUg = rhog*steam.diffProp(steam.steamIntEnergy, p)+Ug*drhog
#    print Uw, Ug, drhoUw, drhoUg
    return [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg];

pi = 600;      # inject pressure: psia
p0 = 500;      # cell pressure: psia
Sg0 = 1.0;     # init saturation
#Hw = 471.6;   # inject water enthalpy: btu/lb
#J = 10.0;     # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");
getcontext().prec = 5;

Hw = 262.09;
Jvalues = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,\
           0.0072, 0.0073, 0.00732, 0.00738, 0.0074, 0.0075, 0.0077,\
           0.008, 0.0085, 0.0087, 0.00875, 0.00879, 0.00881,\
           0.00885, 0.009, 0.01, 0.02, 0.03];
"""
dSg=[]
dp=[]
compressibility=[]
steamSaturation=[]
pressure=[]
residual=[]
"""

p_conv = []; Sg_conv = [];
p1 = []; Sg1 = [];

print '\\begin{table}[H]\n\\centering'
print '\\begin{tabular}{ r r r r r r r r r }'
print 'J & $\\delta S_g$ & $\\delta p$ & $\\hat{a}_{22}$'
print '& $S_g$ & $p$ & $\\hat{R}_e$ & conv $S_g$ & conv $p$ \\\\'

for J in Jvalues:
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
    print '& $S_g$ & $p$ & $\\hat{R}_e$'
    """
    Sg = Sg0;
    p = p0;
    [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
    Rw0 = rhow*(1-Sg) + rhog*Sg;
    Re0 = (rhow*(1-Sg)*Uw + rhog*Sg*Ug);

    for iter in range(1):
        Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
        Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
        RHS = -np.array([Rw, Re]).transpose();
        a11 = rhog-rhow;
        a21 = rhog*Ug-rhow*Uw;
        a12 = (1-Sg)*drhow + Sg*drhog + J;
        a22 = (1-Sg)*drhoUw + Sg*drhoUg + J*Hw;
        comp = a22 - a21/a11*a12;
        Rebar = Re - a21/a11*Rw;
        A = np.array([[a11, a12], [a21, a22]]);
        #    print A
        x = np.linalg.solve(A, RHS);
        Sg = Sg + x[0];
        p = p + x[1];
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
    """
    dSg.append(x[0])
    dp.append(x[1])
    compressibility.append(comp)
    steamSaturation.append(Sg)
    pressure.append(p)
    residual.append(Rebar)
    """
    p1.append(p); Sg1.append(Sg);
    print J, '&',\
              Decimal(x[0]).normalize(),'&',\
              Decimal(x[1]).normalize(),'&',\
              Decimal(comp).normalize(),'&',\
              Decimal(Sg).normalize(),'&',\
              Decimal(p).normalize(),'&',\
              Decimal(Rebar).normalize(),'&'
    for iter in range(19):
        Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
        Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
        RHS = -np.array([Rw, Re]).transpose();
        a11 = rhog-rhow;
        a21 = rhog*Ug-rhow*Uw;
        a12 = (1-Sg)*drhow + Sg*drhog + J;
        a22 = (1-Sg)*drhoUw + Sg*drhoUg + J*Hw;
        comp = a22 - a21/a11*a12;
        Rebar = Re - a21/a11*Rw;
        A = np.array([[a11, a12], [a21, a22]]);
        #    print A
        x = np.linalg.solve(A, RHS);
        Sg = Sg + x[0];
        p = p + x[1];
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
    print '    ', Decimal(Sg).normalize(), '&',\
          Decimal(p).normalize(), '\\\\'
    p_conv.append(p); Sg_conv.append(Sg);

print '\\end{tabular}\n\\end{table}'
print '\\newpage'

# visualization
#plt.plot(Jvalues, dSg)
#plt.plot(Jvalues, dp)
#plt.plot(Jvalues, compressibility)
#plt.plot(Jvalues, steamSaturation)
#plt.plot(Jvalues, pressure)
#plt.plot(Jvalues, residual)
plt.figure(figsize=(10, 16))
plt.subplot(2, 1, 1)
plt.plot(Jvalues, p1, 'o-')
plt.plot(Jvalues, p_conv, '^-')
plt.plot([0, 0.03], [0, 0])
plt.plot([0, 0.03], [500, 500], 'k')
plt.legend(['First iteration $p$','Converged $p$','$p=0$','Initial $p=500$'], fontsize=16, loc='lower right')
#plt.xlabel('Injectivity $J$', fontsize=20)
plt.ylabel('Pressure (psi)', fontsize=20)
plt.ylim((-400, 1000))
plt.tick_params(labelsize=16)
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(Jvalues, Sg1, 'o-')
plt.plot(Jvalues, Sg_conv, '^-')
plt.legend(['First iteration $S_g$','Converged $S_g$'], fontsize=16, loc='lower right')
plt.xlabel('Injectivity $J$', fontsize=20)
plt.ylabel('Steam Saturation', fontsize=20)
plt.ylim((0.0, 1.5))
plt.tick_params(labelsize=16)
plt.grid(True)
plt.show()
