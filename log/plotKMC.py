print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)

runID = 3
dt = 5e-4
frame_w = dt
plot_interval = 20

nr_list = [31, 32]
lambda_list = [1e-2, 1e-2]

fig, ax = mpl.pyplot.subplots()

for i in range(len(nr_list)):
    dirpath = os.getcwd().removesuffix("/ld_analysis/log") + \
                "/LD/LD-cpp/Data{nr}".format(nr = nr_list[i])
    fpath = dirpath + "/kmc.dat"

    [step, numA, numB] = file_utils.readData(fpath, 1)
    step = [i*frame_w for i in step]
    step = step[0:len(step):plot_interval]
    numA = numA[0:len(numA):plot_interval]
    numB = numB[0:len(numB):plot_interval]
    
    y = [numB[j]/numA[j] for j in range(len(step))]

    ax.set(xlim = (0,15e3))
    ax.set_xlabel(r"$t~(\sigma^2/D)$", fontsize = 14)
    ax.set_ylabel(r"$K = N_B/N_A$", fontsize = 14)
    ax.plot(step, y, lw = 1, ls = "-", label = r"$\lambda =$" + str(lambda_list[i]))

ax.axhline(y=1,xmin=0,xmax=1, lw = 0.5, color="k", ls = "-")
ax.axhline(y=0,xmin=0,xmax=1, lw = 0.5, color="k", ls = "-")
ax.legend(fontsize=10)

figpath = "../../LD/LD-cpp/imgs/KMCplot2.png"
#fig.savefig(figpath, dpi = 600, bbox_inches='tight')