import numpy as np
import matplotlib.pyplot as plt
import steamProp

class ConnectionList:
    def __init__(self):
        self._connection = [(1,),(0,)];
        """
        for i in range(1):
            self._connection.append((i - 1, i + 1));
        """

    def GetConnection(self, i):
        return self._connection[i];

    def GeoTrans(self, reservoir, i, j):
        ki = reservoir.Permeability(i);
        kj = reservoir.Permeability(j);
        dxi = reservoir.Deltax(i);
        dxj = reservoir.Deltax(j);
        km = (dxi + dxj)/(dxi/ki + dxj/kj);
        return km * reservoir.GetSectionA() * 2 / (dxi + dxj);

    def FluidTrans(self, reservoir, pressureField, \
                   saturationField, i, j, diffVar=0):
        u = self.Upstream(pressureField, i, j);
        if diffVar and u != i:
            return 0;
        S = [saturationField[u], 1-saturationField[u]];
        p = pressureField[u];
        kr = np.array(S);
        steam = reservoir.getFluid();
        muw = steam.waterViscosity(p);
        mug = steam.steamViscosity(p);
        mu = np.array([mug, muw]);
        if diffVar == 2:
            dmuw = steam.diffProp(steam.waterViscosity, p);
            dmug = steam.diffProp(steam.steamViscosity, p);
            dmu = np.array([dmug, dmuw]);
            return -kr*dmu/np.power(mu, 2);
        return {0:kr/mu, 1:np.array([1,-1])/mu}[diffVar];

    def Transmissibility(self, reservoir, pressureField, \
                         saturationField, i, j, diffVar=0):
        T = self.GeoTrans(reservoir, i, j);
        l = self.FluidTrans(reservoir, pressureField, \
                            saturationField, i, j, diffVar);
        return T*l;

    def HeatTrans(self, reservoir, pressureField, \
                  saturationField, i, j, diffVar=0):
        u = self.Upstream(pressureField, i, j);
        if diffVar and u!=i:
            return 0;
        steam = reservoir.getFluid();
        p = pressureField[u];
        if diffVar == 2:
            dHw = steam.diffProp(steam.waterEnthalpy, p);
            dHg = steam.diffProp(steam.steamEnthalpy, p);
            dH = np.array([dHg, dHw]);
            return dH;
        Hw = steam.waterEnthalpy(p);
        Hg = steam.steamEnthalpy(p);
        H = np.array([Hg, Hw]);
        return {0:H, 1:0}[diffVar];

    def Upstream(self, pressureField, i, j):
        if pressureField[i] > pressureField[j]:
            return i;
        return j;


class Reservoir:
    def __init__(self):
        self._fluid = steamProp.steamProp('saturated_steam.org');
        self._size = 2;
        self._nPrimVar = 2;
        self._poreVol = [1000.0, 1000.0];
        self._perm = [20.9, 20.9];
        self._deltax = 10.0;
        self._sectionA = 100.0;
        self._connList = ConnectionList();

    def Size(self):
        return self._size;

    def GetConnection(self, i):
        return self._connList.GetConnection(i);
    def Deltax(self, i):
        return self._deltax;
    def GetSectionA(self):
        return self._sectionA;

    def CheckBound(self, i):
        return (i >= 0) and (i < self._size);
    def PoreVolumeList(self):
        return self._poreVol;
    def PoreVolume(self, i):
        if (self.CheckBound(i)):
            return self._poreVol[i];
        else:
            print 'Coordinate out of bound'
            return -1;
    def PermeabilityList(self):
        return self._perm;
    def Permeability(self, i):
        if (self.CheckBound(i)):
            return self._perm[i];
        else:
            print 'Coordinate out of bound'
            return -1;
    def getFluid(self):
        return self._fluid;
    def Transmissibility(self, i, j, x, diffVar=0):
        #adding i j order
        n = self.Size();
        Sg = x[0:n];
        p = x[n:2*n];
        T = self._connList.Transmissibility(self, p, Sg, i, j);
        H = self._connList.HeatTrans(self, p, Sg, i, j);
        if diffVar:
            dT = self._connList.Transmissibility(self, p, Sg, i, j, diffVar);
            dH = self._connList.HeatTrans(self, p, Sg, i, j, diffVar);
            return (0, 0); #return (dT, T*dH+H*dT);
        return (1.99e3, 9.38683e5); #return (T, T*H);