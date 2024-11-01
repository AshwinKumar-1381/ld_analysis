print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m

mpl.pyplot.rc('text', usetex=True)

def readLogParam(entry, fpath):
    fobj = open(fpath, mode = "r", encoding = "utf-8")
    for line in fobj:
        line = line.removesuffix("\n").split(sep = " ")
        if(line[0] == entry):
            return(liteval(line[1]))
            break
    fobj.close()

def makeFilePath(nr, type1, type2, extras):
    fpath = os.getcwd().removesuffix("/ld_analysis/log") + \
            "/LD/{type1}/Data{nr}".format(type1 = type1, nr = nr)
    if(extras[0] == 'log'): fpath += "/log.dat"
    elif(extras[0] == 'dump'): 
        fpath += "/relaxa_{type2}/dump_{type2}.dat".format(type2 = type2)
    elif(extras[0] == 'thermo'): fpath += "/thermo.dat"
    elif(extras[0] == 'fig'): fpath = fpath.removesuffix("/{type1}/Data{nr}".format(type1 = type1, nr = nr)) +\
                                      "/imgs/msd_{type1}.png".format(type1=type1)
    elif(extras[0] == 'msd' and type1 == 'cpp'): fpath += "/relaxa_{type2}/msd.dat".format(type2 = type2) 
    return(fpath)

def readData(fpath):
    fobj = open(fpath, mode = 'r', encoding = "utf-8")
    nCols = len(fobj.readline().removesuffix("\n").split(sep = " "))
    Cols = [[] for i in range(nCols)]
    
    for line in fobj:
        line = line.removesuffix("\n").split(sep = " ")
        for i in range(nCols): Cols[i].append(liteval(line[i]))
    
    fobj.close()
    return(Cols)

colors = ['C0', 'C1', 'C2', 'C3', 'C4']
nr2 = [1,2,3,4,5]
damp2 = [0.01, 0.05, 0.1, 0.5, 1]
dt2 = [0.5e-4]*len(nr2)

fig2, ax2 = mpl.pyplot.subplots()
for i in range(len(nr2)):
    
    fpath = makeFilePath(nr2[i], 'lmp', 'ss', ['thermo'])
    
    [t, T, pe, msd] = readData(fpath)
    t = [j*dt2[i] for j in t[1:len(t)-1]]
    msd = [j for j in msd[1:len(msd)-1]]
    
    print((m.log(max(msd))-m.log(min(msd)))/(m.log(max(t))-m.log(min(t))))
    
    ax2.loglog(t, msd, linestyle = '-', lw = 1.2, label = str(damp2[i]))


fig2.legend(title = r'$D^*~(lmp)$', bbox_to_anchor = (0.1,0.8,0.13,0.08))
ax2.set(xlim=(8,1.2e4))
ax2.set_xlabel("$t_2^*~(\sigma\sqrt{m/\epsilon})$", fontsize = 15)
ax2.set_ylabel("MSD", fontsize = 12)

"""
ax2t = ax2.twiny()
nr1 = [6,7,8,9,10]
damp1 = [1e-4, 2.5e-3,0.01,0.25,1]
dt1 = [5e-7, 2.5e-6, 5e-6, 2.5e-5, 5e-5]

for i in range(len(nr1)):
    
    fpath = makeFilePath(nr1[i], 'cpp', 'eq', ['msd'])
    [t, msd] = readData(fpath)

    print(len(t))
    t = [j/m.sqrt(damp1[i]) for j in t[1000:len(t)-1]]
    msd = [j for j in msd[1000:len(msd)-1]]
    
    ax2t.loglog(t, msd, linestyle = '--', color = colors[i], lw = 0.9, label = str(damp1[i]))

ax2t.legend(title = r"$m^*~(cpp)$", bbox_to_anchor = (0.1,0.8,0.18,0.2))


ax2t.set_xlabel("$t_2^*~(\sigma\sqrt{m/\epsilon})$", fontsize = 15)
ax2t.set(xlim=(8,1.2e4))
"""

fpath = makeFilePath(6, 'lmp','eq',['fig'])
#fig2.savefig(fpath, dpi = 600, bbox_inches='tight')
