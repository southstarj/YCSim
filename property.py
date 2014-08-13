import numpy as np
import scipy as sp
import matplotlib as plt

"""
p = np.zeros(0);
bp = np.zeros(0);
rho = np.zeros(0);

with open('prop_steam.txt', 'r') as f:
    line = f.readline();
    f.readline();
    while True:
        line = f.readline();
        values = [float(x) for x in line.split()];
        if len(values) == 0:
            break;
        tD_ref = np.append(tD_ref, values[3]*1e3/(Dx*Dy*Dz*phi/qT));
        qt_ref = np.append(qt_ref, values[2]);
        qp_ref = np.append(qp_ref, values[4]+values[5]);
        Fw_ref = np.append(Fw_ref, values[5]/values[2]);
        Fp_ref = np.append(Fp_ref, values[5]/(values[4]+values[5]));
        Pm_ref = np.append(Pm_ref, values[9]);
"""
