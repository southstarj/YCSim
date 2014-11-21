import numpy as np
import matplotlib.pyplot as plt
import steamProp

class ConnectionList:
    def __init__(self):
        self._connection = [()];
        for i in range(1, 9):
            self._connection.append((i - 1, i + 1));
        self._connection.append(());

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
        if diffVar and u!=i:
            return 0;
        S = [saturationField[u], 1-saturationField[u]];
        p = pressureField[u];
        kr = np.array(S);
        steam = reservoir.getFluid();
        if diffVar == 2:
            dmuw = steam.diffProp(steam.waterViscosity, p);
            dmug = steam.diffProp(steam.steamViscosity, p);
            dmu = np.array([dmug, dmuw]);
        muw = steam.waterViscosity(p);
        mug = steam.steamViscosity(p);
        mu = np.array([mug, muw]);
        return {0:kr/mu, 1:1/mu, 2:-kr*dmu/np.power(mu, 2)}[diffVar];

    def Transmissibility(self, reservoir, pressureField, \
                         saturationField, i, j, diffVar=0):
        T = GeoTrans(reservoir, i, j);
        l = FluidTrans(reservoir, pressureField, \
                       saturationField, i, j, diffVar);
        return T*l;

    def HeatTrans(self, reservoir, pressureField, \
                  saturationField, i, j, diffVar=0):
        u = self.Upstream(pressureField, i, j);
        if diffVar and u!=i:
            return 0;
        steam = reservoir.getFluid();
        if diffVar == 2:
            dHw = steam.diffProp(steam.waterEnthalpy, p);
            dHg = steam.diffProp(steam.steamEnthalpy, p);
            dH = np.array([dHg, dHw]);
        Hw = steam.waterEnthalpy(p);
        Hg = steam.steamEnthalpy(p);
        H = np.array([Hg, Hw]);
        return {0:H, 1:0, 2:dH}[diffVar];

    def Upstream(self, pressureField, i, j):
        if pressureField[i] > pressureField[j]:
            return i;
        return j;


class Reservoir:
    def __init__(self):
        self._fluid = steamProp.steamProp('saturated_steam.org');
        self._size = 10;
        self._nPrimVar = 2;
        self._poreVol = [1000 for i in range(10)];
        self._perm = [20.9 for i in range(10)];
        self._deltax = 10;
        self._sectionA = 100;
        self._connList = ConnectionList();

    def Size(self):
        return self._size;

    def GetConnection(self, i):
        return self._connList.GetConnection(i);

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
        T = _connList.Transmissibility(self, p, Sg, i, j);
        H = _connList.HeatTrans(self, p, Sg, i, j);
        if diffVar:
            dT = _connList.Transmissibility(self, p, Sg, i, j, diffVar);
            dH = _connList.HeatTrans(self, p, Sg, i, j, diffVar);
            return (dT, T*dH+H*dT);
        return (T, T*H);