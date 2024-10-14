"""
 Code to read time-averaged Particle Density Function (PDF)
 data corresponding to multiple runs
 Author: Ashwin Kumar M (CH23S006)
 Date created: 09.09.24
 Last modified: 09.09.24
"""

print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
from ast import literal_eval as liteval
import os

mpl.pyplot.rc('text', usetex=True)

nr_list = [3, 1, 4, 5]
dir = ["-", "DOWN", "UP", "UP"]
markers = ["v", '^', "s", ">"]
line_style = ["-", "-", "-", "--"]

lw = 0.4
cblue = "#2D51F5"
cred = "#F34943"

fig, ax = mpl.pyplot.subplots()

x = []
vel = []
for i in range(len(nr_list)):
    fpath = os.getcwd().removesuffix("/code") + \
            "/LD/Data{nr}/relaxa_ss/current_{nr}.dat".format(nr = nr_list[i])
    
    f_obj = open(fpath, mode = "r", encoding = "utf-8")
    
    line = f_obj.readline().split(sep = " ")
    x_start = liteval(line[0])
    x_end = liteval(line[1])
    n_bin = liteval(line[2])
    L = liteval(line[3])
    N = liteval(line[4])
    PeA = liteval(line[5])
    PeB = liteval(line[6].removesuffix("\n"))
    
    for j in range(n_bin):
        line = f_obj.readline().split(sep = " ")
        x.append(x_start + liteval(line[0])*(x_end - x_start)/n_bin)
        vel.append(liteval(line[1].removesuffix("\n")))
    
    ax.plot(x, vel, ls = line_style[i], lw = lw, color = "k", marker = markers[i], ms = 3, 
            mec = "k", mfc = "white", label = "$Pe_A = {PeA}, Pe_B = {PeB}, {dir}$".
            format(PeA = PeA, PeB = PeB, dir = dir[i]))
    ax.set_xlim(x_start, x_end + 1)
    ax.set_ylim(-1e-3, 1e-3)
    ax.set_xlabel("x ($\sigma$)")
    ax.set_ylabel(r"Particle current ($\sigma / \tau$)")
    ax.legend(fontsize = 6)

    x = []
    vel = []
    
    f_obj.close()

#fig.savefig("particle_current.png", dpi = 1000, bbox_inches = 'tight')