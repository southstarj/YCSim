import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

class steamProp:
    def __init__(self, filename):
        ncol = 10;
        self.data = np.zeros((40, ncol));
        # read steam propertie table from file
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
                i += 1;
        # set up spline interpolation object
        self.tck = [];
        for i in range(ncol):
            self.tck.append(interpolate.splrep(self.data[:,0],\
                                               self.data[:,i], s=0));

    def extractData(self, p, colNum):
    # Get data from interpolation
        if type(p) is list:
            _value = [];
            for i in range(len(p)):
                _value[i] = interpolate.splev(p, self.tck[colNum], der=0);
        else:
            _value = interpolate.splev(p, self.tck[colNum], der=0);
        return _value;

    def diffProp(self, propFunc, p):
    # Numerical property function differentiation
        dp = 1;
        if type(p) is list:
            _value = [];
            for i in range(len(p)):
                rightDp = propFunc(p[i] + dp);
                leftDp = propFunc(p[i] - dp);
                _value[i] = (rightDp - leftDp) / (2 * dp);
        else:
            rightDp = propFunc(p + dp);
            leftDp = propFunc(p - dp);
            _value = (rightDp - leftDp) / (2 * dp);
        return _value;

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
    def waterViscosity(self, p):
        return self.extractData(p, 8);
    def steamViscosity(self, p):
        return self.extractData(p, 9);

