"""
 Code to read time-averaged Particle Density Function (PDF)
 data corresponding to multiple runs
 Author         : Ashwin Kumar M (CH23S006)
 Date created   : 09.09.24
 Last modified  : 16.10.24
"""

print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
from ast import literal_eval as liteval
import os

mpl.pyplot.rc('text', usetex=True)

def makeFilePath(nr, type1, type2, extras):
    fpath = os.getcwd().removesuffix("/ld_analysis/log") + \
            "/LD/{type1}/Data{nr}".format(type1 = type1, nr = nr)
            
    if(extras[0] == 'current'):
        if(type1 == 'lmp'): fpath += "/current{nr}.dat".format(nr = nr)
        elif(type1 == 'cpp'): 
            fpath += "/relaxa_{type2}/current{nr}.dat".format(type2 = type2, nr = nr)
            
    elif(extras[0] == 'fig'):
        fpath = fpath.removesuffix("/Data{nr}".format(nr = nr))
        fpath += "/imgs/currents10.png"
    return(fpath)

nr = [7, 9, 10, 11, 8]
nr = [10]
type1 = 'lmp'
type2 = 'ss'

lw = 0.6

fig, ax = mpl.pyplot.subplots()
ax.axhline(y = 0, xmin = 0, xmax = 1, color = 'k', linestyle = '-', lw = 0.4) 
ax.axvline(x = 35, ymin = 0, ymax = 1, color = 'k', linestyle = '-', lw = 0.4, alpha = 0.8)
ax.axvline(x = 65, ymin = 0, ymax = 1, color = 'k', linestyle = '-', lw = 0.4, alpha = 0.8)

x = []
vel = []
for i in range(len(nr)):
    fpath = makeFilePath(nr[i], type1, type2, ['current'])
    
    fobj = open(fpath, mode = "r", encoding = "utf-8")
    
    line = fobj.readline().removesuffix("\n").split(sep = " ")
    line = [liteval(i) for i in line]
    [x_start, x_end, n_bin, L, N, PeA, PeB] = line
    
    line = fobj.readline()
    
    for j in range(n_bin):
        line = fobj.readline().removesuffix("\n").split(sep = " ")
        x.append(x_start + liteval(line[0])*(x_end - x_start)/n_bin - 0.5)
        vel.append(liteval(line[1]))
    
    ax.plot(x, vel, lw = lw, marker = "o", ms = 3, 
            label = "$Pe_A = {PeA}, Pe_B = {PeB}$".format(PeA = PeA, PeB = PeB))
    ax.set_xlim(x_start, x_end)
    ax.set_ylim(-2e-2, 2e-2)
    ax.set_xlabel("x ($\sigma$)", fontsize = 13)
    ax.set_ylabel(r"Particle current ($\sigma / \tau$)", fontsize = 13)
    ax.set_xticks(ticks = range(int(x_start),int(x_end+2),10))
    ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend(fontsize = 10)

    x = []
    vel = []
    
    fobj.close()

figPath = makeFilePath(0, type1, type2, ['fig'])
#fig.savefig(figPath, dpi = 600, bbox_inches = 'tight')