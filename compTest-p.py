import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import steamProp as prop
import Reservoir
import Fluid
from Simulator import *

n = 1;
nv = 2;
debugstep = -1;
# geological settings
poreVol = [1000.0 for i in range(n)];
perm = [0.5 for i in range(n)];
deltax = 100.0;
Area = 10.0;
reservoir = Reservoir.Reservoir(n, poreVol, perm, deltax, Area);
# initial distribution
p0 = [500 for i in range(n)];  p0[0] = 300;
Sg0 = [1.0 for i in range(n)]; Sg0[0] = 1.0;
Wellnum = 0; qT = -5.0; pWater = 1; pInj = 300;          # water injection

rhoB = reservoir.getFluid().waterDensity(pInj);
HB = reservoir.getFluid().waterEnthalpy(pWater);
muB = reservoir.getFluid().waterViscosity(pWater);
T = reservoir.Permeability(Wellnum)*reservoir.GetSectionA()\
    /(reservoir.Deltax(Wellnum)*muB);

dt = 1;                           # time step

np.set_printoptions(precision=5);

def compressibility_Test(reservoir, p, Sg, Hw):
    steam = reservoir.getFluid();
    rhog = steam.steamDensity(p);
    rhow = steam.waterDensity(p);
    Ug = steam.steamIntEnergy(p);
    Uw = steam.waterIntEnergy(p);
    #print '  rho, U =', rhog, rhow, Ug, Uw;
    beta = (rhog*Ug - rhow*Uw)/(rhog-rhow);
    #print '  beta =', beta;
    drhow = steam.diffProp(steam.waterDensity, p);
    drhog = steam.diffProp(steam.steamDensity, p);
    dUw = steam.diffProp(steam.waterIntEnergy, p);
    dUg = steam.diffProp(steam.steamIntEnergy, p);
    drhoUw = np.array([(x*dy+y*dx) for x,y,dx,dy in zip(rhow,Uw,drhow,dUw)]);
    drhoUg = np.array([(x*dy+y*dx) for x,y,dx,dy in zip(rhog,Ug,drhog,dUg)]);
    Sw = 1-Sg;
    alpha1 = (beta - Hw);
    alpha2 = (beta*(Sw*drhow+Sg*drhog)-(Sw*drhoUw+Sg*drhoUg));
    print '  alpha =', alpha1/alpha2;
    return alpha1/alpha2, beta;

#print 'Hw =', HB
Sg0 = 1.0
ps = np.array(range(1,299))#[60, 100, 150, 200, 250, 290, 295, 300]);
pwaters = np.array([pWater]);
HB = reservoir.getFluid().waterEnthalpy(pwaters);

alpha, beta = compressibility_Test(reservoir, ps, Sg0, HB);

lincomp = np.array([]); dp = np.array([]);
#p = 300;
for p in ps:
    x = np.array([Sg0, p]);
    x0 = x;
    RHS = GenerateRHS(reservoir, dt, x, x0);
    #RHS = BoundaryCond_Rate(reservoir, RHS, qT, pWater, pInj, Wellnum);
    A = GenerateJacobian(reservoir, dt, x);
    A, RHS = BoundaryCond_Pres(reservoir, A, RHS, pwaters, pInj,\
             x[reservoir.Size()+Wellnum], Wellnum);
    dx, a22bar, Rebar = LinearSolver(reservoir, dt, A, RHS);
    
    lincomp = np.append(lincomp, -Rebar/(a22bar*T*(pInj-p)));
    dp = np.append(dp, p-Rebar/(a22bar));
    
    """
    print 'p =', p
    compressibility_Test(reservoir, p, 0.1, HB);
    print '  Injectivity =', T*(pInj - p)/poreVol[0]
    print '  A ='
    print A
    print '  RHS =', RHS
    print '  dx =', dx[0], dx[1]
    print '  comp =', comp
    print '  Rebar =', Rebar
    """

#print 'a22bar =', a22bar
fig = plt.figure()#figsize=(16, 9), dpi = 80);
Axbeta = fig.add_subplot(311)
plt.title('$S_g = '+str(Sg0)+', \; H_w^0 = 28.08Btu/lbm(60F^o)$')
Axalpha = fig.add_subplot(312)#, sharex=Axbeta)
Axcomp = fig.add_subplot(313, sharex=Axalpha)
p1, = Axbeta.plot(ps, beta)
Axbeta.grid(True)
Axbeta.legend([p1], ['$\\beta$'], loc=4)
p1, = Axalpha.plot(ps, alpha)
p2, = Axalpha.plot(ps, dp)
Axalpha.grid(True)
#Axalpha.set_ylim((-300, 300))
Axalpha.legend([p1,p2],['$\\alpha$', '$p+\\delta p$'], loc=7)
p1, = Axcomp.plot(ps, lincomp)
Axcomp.grid(True)
Axcomp.legend([p1],['$\\hat{\\alpha}$'])
plt.show()
