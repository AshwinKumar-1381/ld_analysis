print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
from file_utils import *

mpl.pyplot.rc('text', usetex=True)

nr2 = 2
fpath2 = makeFilePath(nr2,'lmp','eq',['msd'])
[t2, msd2] = readData(fpath2)
t2 = t2[1:len(t2)]
msd2 = msd2[1:len(msd2)]

avg_dev = []
nr1 = [7, 19, 20, 21, 22, 23, 24, 25, 26]
for nr in nr1:
    del_msd = []
    
    fpath1 = makeFilePath(nr, 'cpp', 'eq', ['msd'])
    [t1, msd1] = readData(fpath1)
    t1 = t1[1:len(t1)]
    msd1 = msd1[1:len(msd1)]
    
    for i in range(len(t1)):
        del_msd.append(abs(msd2[i]-msd1[i])/msd2[i]*100)
    
    avg_dev.append(round(sum(del_msd)/len(del_msd),3))
print(avg_dev)

timestep = [2.5e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
timestep = [str(i) for i in timestep]

fig, ax = mpl.pyplot.subplots()
ax.semilogy(timestep, avg_dev, 'o-', lw = 0.8, ms = 4)
ax.set_xlabel(r"cpp timestep")
ax.set_ylabel(r"$<\delta MSD>_t (\%)$")

#fig.savefig("del_msd.png", dpi = 600, bbox_inches = 'tight')