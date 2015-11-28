import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
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

def calcComp(p, Sg):
    [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
    

pi = 600;      # inject pressure: psia
p0 = 500;      # cell pressure: psia
Sg0 = 1.0;     # init saturation
Hw = 262.25;   # inject water enthalpy: btu/lb
#J = 0.003;      # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");

#Hwvalues = range(400);
#pvalues = range(299);
Jvalues = np.linspace(0.0001, 0.014, 1000);
alphaList = np.array([]);
comp_list = np.array([]);
dp_list = np.array([]);

for J in Jvalues:
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
        Sjdrhoj = (1-Sg)*drhow + Sg*drhog;
        SjdrhoUj =(1-Sg)*drhoUw + Sg*drhoUg;
        a12 = Sjdrhoj + J;
        a22 = SjdrhoUj + J*Hw;
        beta = a21/a11;
        alpha1 = beta - Hw;
        alpha2 = beta*Sjdrhoj - SjdrhoUj;
        alpha = alpha1/alpha2;
        alphaList = np.append(alphaList, alpha);

        a22bar = a22 - beta*a12;
        Rebar = Re - beta*Rw;
        delta_p = -Rebar/a22bar;
        comp_linear = delta_p/(J*(pi - p));
        comp_list = np.append(comp_list, comp_linear);

        A = np.array([[a11, a12], [a21, a22]]);
        #    print A
        x = np.linalg.solve(A, RHS);
        Sg = Sg + x[0];
        p = p + x[1];
        dp_list = np.append(dp_list, p);
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);

    #print '    ', Sg, p

plt.figure(figsize=(16, 12))
plt.plot(Jvalues, alphaList)
plt.plot(Jvalues, comp_list)
plt.plot(Jvalues, dp_list)
plt.xlabel('Injectivity $J$', fontsize=25)
plt.tick_params(labelsize=20)
plt.ylim(-1000,1000)
plt.grid(True)
plt.legend(['Nonlinear compressibility', 'Linear compressibility', 'Pressure after first iteration'],loc='lower left',fontsize=20)
plt.show()
