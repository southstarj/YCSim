import numpy as np
import matplotlib.pyplot as plt
import steamProp

class ConnectionList:
    def __init__(self):
        self._connection = [()];
        for i in range(1, 9):
            self._connection.append((i - 1, i + 1));
        self._connection.append(());
        
    def GeoTrans(self, reservoir, i, j):
        ki = reservoir.Permeability(i);
        kj = reservoir.Permeability(j);
        dxi = reservoir.Deltax(i);
        dxj = reservoir.Deltax(j);
        km = (dxi + dxj)/(dxi/ki + dxj/kj);
        return km * reservoir.GetSectionA() * 2 / (dxi + dxj);

    def FluidTrans(self, reservoir, i, j):
        


class Reservoir:
    def __init__(self):
        self._fluid = steamProp.steamProp('saturated_steam.org');
        self._size = 10;
        self._nPrimVar = 2;
        self._poreVol = [1000 for i in range(10)];
        self._perm = [20.9 for i in range(10)];
        self._deltax = 10;
        self._sectionA = 100;

    def Size(self):
        return self._size;
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
    def Transmissibility(self, i, p, Sg):
        k = self.Permeability(i);
        steam = self._fluid;
        muw = steam.waterViscosity(p);
        mug = steam.steamViscosity(p);
        Tw = k/muw; Tg = k/mug;
        return (Tw*(1-Sg), Tg*Sg);