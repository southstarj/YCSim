import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

class steamProp:
    def __init__(self, filename):
        self.data = np.zeros((40, 8));
        with open(filename, 'r') as f:
            line = f.readline();
            f.readline();
            i = 0;
            while True:
                line = f.readline();
                values = [];
                for x in line.split():
                    if x == '|':
                        continue;
                    values.append(float(x));
                if len(values) == 0:
                    break;
                self.data[i] = values;
                i = i + 1;
        self.tck = [];
        for i in range(8):
            self.tck.append(interpolate.splrep(self.data[:,0], self.data[:,i], s=0));

    def extractData(self, p, colNum):
        return interpolate.splev(p, self.tck[colNum], der=0);
    def diffProp(self, propFunc, p):
        dp = 1;
        rightDp = propFunc(p + dp);
        leftDp = propFunc(p - dp);
        return (rightDp - leftDp) / (2 * dp);

    def boilingPoint(self, p):
        return self.extractData(p, 1);
    def waterDensity(self, p):
        return 1.0/self.extractData(p, 2);
    def steamDensity(self, p):
        return 1.0/self.extractData(p, 3);
    def waterEnthalpy(self, p):
        return self.extractData(p, 4);
    def steamEnthalpy(self, p):
        return self.extractData(p, 5);
    def waterIntEnergy(self, p):
        return self.extractData(p, 6);
    def steamIntEnergy(self, p):
        return self.extractData(p, 7);
