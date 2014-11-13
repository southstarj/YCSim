import numpy as np
import matplotlib.pyplot as plt
import steamProp

class Reservoir:
    def __init__(self, n, Vp, k):
        self._fluid = steamProp.steamProp('saturated_steam.org');
        self._size = n;
        self._poreVol = Vp;
        self._perm = k;

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