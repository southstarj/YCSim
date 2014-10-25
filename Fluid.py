"""
class Fluid: only for water/steam now
created Oct.24 2014 by Yifan Wang
"""

class Fluid:
    def __init__(self):
        self._steam = prop.steamProp('saturated_steam.org')
        self._pressure = 0;

    def computeViscosity(self):
        return 0
    def calcProp(self, p):
        if (p == self._pressure):
            return
        steam = self._steam
        crho = 0.062428*0.178107607;    # kg/m3 to lb/RB
        self._rho = [steam.waterDensity(p),\
                     steam.steamDensity(p)]
        self._drho = [steam.diffProp(steam.waterDensity, p),\
                      steam.diffProp(steam.steamDensity, p)]
        cU = 0.429923;      # kJ/kg to Btu/lb
        self._U = [steam.waterIntEnergy,\
                   steam.steamIntEnergy]
        self._drhoU = [rho*steam.diffProp(steam.waterIntEnergy, p)\
                          +Uw*drhow,\
                       rhog*steam.diffProp(steam.steamIntEnergy, p)\
                          +Ug*drhog]

    def computeDensity(self, phaseNum, p):
        self.calcProp(p)
        return self._rho[phaseNum]

    def computeDensityDiff(self, phaseNum, p):
        self.calcProp(p)
        return self._drho[phaseNum]

    def computeInternalEnergy(self, phaseNum, p):
        self.calcProp(p)
        return self._U[phaseNum]

    def computeInternalEnergyDiff(self, phaseNum, p):
        self.calcProp(p)
        return self._drhoU[phaseNum]
