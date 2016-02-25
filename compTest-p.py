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
#Hw = 471.6;   # inject water enthalpy: btu/lb
J = 0.1;      # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");

Hw = 262.09;
pvalues = np.linspace(0, 299, 2000);
alphaList = np.array([]);
comp_list = np.array([]);
dp_list = np.array([]);
pconv_list = np.array([]);

for i_p in pvalues:
    Sg = 1.0;
    p = i_p;
    [rhow0, rhog0, drhow0, drhog0, Uw0, Ug0, drhoUw0, drhoUg0] = calcProp(steam, p0);
    Rw0 = rhow0*(1-Sg0) + rhog0*Sg0;
    Re0 = (rhow0*(1-Sg0)*Uw0 + rhog0*Sg0*Ug0);
    [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);

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
        
    for iter in range(10):
        Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
        Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
        RHS = -np.array([Rw, Re]).transpose();
        a11 = rhog-rhow;
        a21 = rhog*Ug-rhow*Uw;
        Sjdrhoj = (1-Sg)*drhow + Sg*drhog;
        SjdrhoUj =(1-Sg)*drhoUw + Sg*drhoUg;
        a12 = Sjdrhoj + J;
        a22 = SjdrhoUj + J*Hw;
        #beta = a21/a11;
        #alpha1 = beta - Hw;
        #alpha2 = beta*Sjdrhoj - SjdrhoUj;
        #alpha = alpha1/alpha2;
        #alphaList = np.append(alphaList, alpha);

        #a22bar = a22 - beta*a12;
        #Rebar = Re - beta*Rw;
        #delta_p = -Rebar/a22bar;
        #comp_linear = delta_p/(J*(pi - p));
        #comp_list = np.append(comp_list, comp_linear);

        A = np.array([[a11, a12], [a21, a22]]);
        #    print A
        x = np.linalg.solve(A, RHS);
        Sg = Sg + x[0];
        p = p + x[1];
        #dp_list = np.append(dp_list, p);
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);

    pconv_list = np.append(pconv_list, p);
    #print '    ', Sg, p

plt.figure(figsize=(16,12))
plt.subplot(211)
plt.plot(pvalues, alphaList, linewidth=2)
plt.plot(pvalues, comp_list, '--', linewidth = 2)
plt.plot([0, 600], [0,0], 'k')
plt.ylim(-200, 200)
plt.xlim(0,300)
plt.tick_params(labelsize=18)
plt.legend(['$\\alpha$', '$\\hat{\\alpha}$'], fontsize=18, loc='lower right')
plt.subplot(212)
plt.plot(pvalues, dp_list, '--', linewidth=2)
plt.plot(pvalues, pconv_list, linewidth=2)
plt.plot([0, 600], [0,0], 'k')
plt.ylim(-200, 800)
plt.xlim(0,300)
plt.tick_params(labelsize=18)
#plt.xlabel('Pressure initial guess(psi)', fontsize=25)
plt.legend(['$p^{k=1}$', '$p_{converge}$'], fontsize=18, loc='lower right')
#plt.title('Injectivity $J=0.1$; Steam saturation initial guess $S_g=0.4$', fontsize=25)
#plt.grid(True)
plt.show()
