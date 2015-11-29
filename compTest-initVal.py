import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
J = 0.1;      # normalized injectivity: lb/cf.psi
steam = prop.steamProp("saturated_steam.org");

pvalues = np.arange(1, 700, 10);
Sgvalues = np.arange(0.0, 1.2, 0.1);
Y, X = np.meshgrid(Sgvalues, pvalues);
print Y.shape, X.shape
alphaList = np.zeros((len(pvalues), len(Sgvalues)));
comp_list = np.zeros((len(pvalues), len(Sgvalues)));
dp_list = np.zeros((len(pvalues), len(Sgvalues)));
res_list = np.zeros((len(pvalues), len(Sgvalues)));
print res_list.shape

for i_Sg in range(len(Sgvalues)):
    for i_p in range(len(pvalues)):
        Sg = Sgvalues[i_Sg];
        p = pvalues[i_p];
        [rhow0, rhog0, drhow0, drhog0, Uw0, Ug0, drhoUw0, drhoUg0] = calcProp(steam, p0);
        Rw0 = rhow0*(1-Sg0) + rhog0*Sg0;
        Re0 = (rhow0*(1-Sg0)*Uw0 + rhog0*Sg0*Ug0);
        [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);

        for iter in range(1):
            Rw = (rhow*(1-Sg) + rhog*Sg) - Rw0 - J*(pi-p);
            Re = (rhow*(1-Sg)*Uw + rhog*Sg*Ug) - Re0 - J*Hw*(pi-p);
            res_list[i_p][i_Sg] = np.sqrt(Rw**2+Re**2);
            
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
            alphaList[i_p][i_Sg] = alpha1/alpha2;

            a22bar = a22 - beta*a12;
            Rebar = Re - beta*Rw;
            delta_p = -Rebar/a22bar;
            comp_linear = delta_p/(J*(pi - p));
            comp_list[i_p][i_Sg] = comp_linear;
            """
            A = np.array([[a11, a12], [a21, a22]]);
            #    print A
            x = np.linalg.solve(A, RHS);
            Sg = Sg + x[0];
            p = p + x[1];
            dp_list = np.append(dp_list, p);
            [rhow, rhog, drhow, drhog, Uw, Ug, drhoUw, drhoUg] = calcProp(steam, p);
            """

        #print '    ', Sg, p

#plt.plot(pvalues, alphaList)
#plt.plot(pvalues, comp_list)
#plt.plot(pvalues, dp_list)
#plt.plot(pvalues, res_list)
#plt.ylim(-300, 300)
plt.figure(figsize=(16,16))
#matplotlib.rcParams['contour.negative_linestyle']='solid';
#res_list=np.log10(res_list/np.amax(res_list));
CS=plt.contour(X,Y,res_list, [-0.1, -0.3, -0.6, -0.9, -1.2, -1.5, -2.5], colors='k')
#CS=plt.contour(X,Y,comp_list)
#CS=plt.contour(X,Y,alphaList,100)
plt.clabel(CS, inline=1, fontsize=14)
#im = plt.imshow(res_list, interpolation='bilinear', origin='lower', cmap=cm.gray)
plt.tick_params(labelsize=18)
plt.xlabel('Pressure initial guess(psi)', fontsize=25)
plt.ylabel('Steam saturation initial guess', fontsize=25)
plt.grid(True)
plt.show()
