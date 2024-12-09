print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "*", "D"]
lw = 0.5
ms = 3
fs = "none"
fontsize = 14

fpath = "./data/coexistPhi.dat"
nDataSets = 4
lamb_val = [r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$1$"]
phiData = file_utils.readData(fpath, nDataSets)

fig, ax = mpl.pyplot.subplots()

for i in range(nDataSets):
    
    Pe = phiData[3*i + 0]
    phiG = phiData[3*i + 1]
    phiL = phiData[3*i + 2]
    
    ax.semilogy(phiG, Pe, marker = markers[i], lw = lw, ms = ms, label = lamb_val[i])
    ax.semilogy(phiL, Pe, marker = markers[i], lw = lw, ms = ms)

ax.set(xlim = (0,1), ylim = (0.1, 100), aspect = 0.3)
ax.set_xlabel(r"$\phi_G , \phi_L$", fontsize = fontsize)
ax.set_ylabel(r"$Pe_s$", fontsize = fontsize)
ax.legend(fontsize = fontsize - 4)