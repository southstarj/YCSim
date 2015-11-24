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
Hw = 471.6;   # inject water enthalpy: btu/lb
J = 10.0;     # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");
getcontext().prec = 5;

Hw = 262.09;
Jvalues = [0.00001, 0.005, 0.007, 0.008, 0.01];

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
for J in Jvalues:
    
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
    

    Sg = Sg0;
    p = p0;
    [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
    Rw0 = rhow*(1-Sg) + rhog*Sg;
    Re0 = (rhow*(1-Sg)*Uw + rhog*Sg*Ug);
    p_vis = [p0,]; Sg_vis = [Sg0,];
    """
    for iter in range(10):
        Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
        Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
        cw = (1-Sg)*drhow + Sg*drhog
        ce = (1-Sg)*drhoUw + Sg*drhoUg
        #Rw = cw*(p-p0)-J*(pi-p)
        #Re = ce*(p-p0)-J*Hw*(pi-p)
        RHS = -np.array([Rw, Re]).transpose();
        a11 = rhog-rhow;
        a21 = rhog*Ug-rhow*Uw;
        a12 = cw + J #(1-Sg)*drhow + Sg*drhog + J;
        a22 = ce + J*Hw #(1-Sg)*drhoUw + Sg*drhoUg + J*Hw;
        comp = a22 - a21/a11*a12;
        Rebar = Re - a21/a11*Rw;
        A = np.array([[a11, a12], [a21, a22]]);
        #    print A
        x = np.zeros(2) #np.linalg.solve(A, RHS);
        x[1] = -Rebar/comp;
        x[0] = (-Rw-a12*x[1])/a11;
        Sg = Sg + x[0];
        p = p + x[1];
        #[rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
        print iter, '&',\
              Decimal(x[0]).normalize(),'&',\
              Decimal(x[1]).normalize(),'&',\
              Decimal(comp).normalize(),'&',\
              Decimal(Sg).normalize(),'&',\
              Decimal(p).normalize(),'&',\
              Decimal(Rebar).normalize(),'\\\\'
    """
    for iter in range(10):
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
        p_vis.append(p);
        Sg_vis.append(Sg);
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
        print '    ', Decimal(Sg).normalize(), '&',\
                      Decimal(p).normalize(), '\\\\'
    

    print '\\end{tabular}\n\\end{table}'
    print '\\newpage'
    plt.figure(figsize=(16,6))
    plt.subplot(1, 2, 1)
    plt.plot(range(11), p_vis)
    plt.xlabel('Iteration number', fontsize=20)
    plt.ylabel('Block pressure(psi)', fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)
    plt.subplot(1, 2, 2)
    plt.plot(range(11), Sg_vis)
    plt.xlabel('Iteration number', fontsize=20)
    plt.ylabel('Steam saturation', fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)

# visualization
#plt.plot(Jvalues, dSg)
#plt.plot(Jvalues, dp)
#plt.plot(Jvalues, compressibility)
#plt.plot(Jvalues, steamSaturation)
#plt.plot(Jvalues, pressure)
#plt.plot(Jvalues, residual)
plt.show()
