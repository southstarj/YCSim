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
T = rhoB*reservoir.Permeability(Wellnum)*reservoir.GetSectionA()\
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
    #print '  alpha =', alpha1, alpha2;
    return alpha1/alpha2, beta, alpha1, alpha2;

#print 'Hw =', HB
Sg0 = 1.0
ps = np.array([60])#[60, 100, 150, 200, 250, 290, 295, 300]);
pwaters = np.array([400]);
HB = reservoir.getFluid().waterEnthalpy(pwaters);

alpha, beta, alpha1, alpha2 = compressibility_Test(reservoir, ps, Sg0, HB);
#print alpha, beta

lincomp = np.array([]); dp = np.array([]);
a22bar_list = np.array([]); Rebar_list = np.array([]);
p = ps[0];
for pw in pwaters:
    x = np.array([Sg0, p]);
    x0 = x;
    RHS = GenerateRHS(reservoir, dt, x, x0);
    #RHS = BoundaryCond_Rate(reservoir, RHS, qT, pWater, pInj, Wellnum);
    A = GenerateJacobian(reservoir, dt, x);
    A, RHS = BoundaryCond_Pres(reservoir, A, RHS, pw, pInj,\
             x[reservoir.Size()+Wellnum], Wellnum);
    dx, a22bar, Rebar = LinearSolver(reservoir, dt, A, RHS);
    a22bar_list = np.append(a22bar_list, a22bar);
    Rebar_list = np.append(Rebar_list, Rebar);

    lincomp = np.append(lincomp, -Rebar/(a22bar*T*(pInj-p)));
    dp = np.append(dp, p-Rebar/(a22bar));
    #print 'lincomp =', lincomp


#print 'a22bar =', a22bar
fig = plt.figure()#figsize=(16, 9), dpi = 80);
Axdp = fig.add_subplot(311)
plt.title('$S_g = '+str(Sg0)+'\\;p = '+str(p)+'psi\\;J='+str(T)+'\\;\\beta='+str(beta)+'$')
Axalpha = fig.add_subplot(312, sharex=Axdp)
Axcomp = fig.add_subplot(313, sharex=Axalpha)
p1, = Axdp.plot(HB, dp)
Axdp.grid(True)
#Axdp.legend([p1],['$p+\\delta p$'], loc=4)
p1, = Axalpha.plot(HB, alpha)
#print alpha1, alpha2
p2, = Axalpha.plot(HB, alpha1)
#p3, = Axalpha.plot(HB, alpha2)
Axalpha.grid(True)
#Axalpha.set_ylim((-300, 300))
#print dp
#Axalpha.legend([p1],['$\\alpha$'], loc=4)
#p1, = Axcomp.plot(HB, a22bar_list)
#p2, = Axcomp.plot(HB, Rebar_list)
#p1, = Axcomp.plot(HB, -Rebar_list/a22bar_list)
p1, = Axcomp.plot(HB, lincomp)
p2, = Axcomp.plot(HB, alpha/(T*alpha+1))
#Axcomp.set_ylim((-700, 100))
Axcomp.grid(True)
#Axcomp.legend([p1],['$\\hat{\\alpha}$'], loc=4)
plt.show()
