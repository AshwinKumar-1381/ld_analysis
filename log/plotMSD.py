print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval

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
    elif(extras[0] == 'fig'): fpath += "/PEplot{nr}.png".format(nr = nr)
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

nr1 = [1,2,3,4,5]
label1 = [0.01, 0.05, 0.1, 0.5, 1]
dt1 = [0.5e-4]*len(nr1)

fig1, ax1 = mpl.pyplot.subplots()
for i in range(len(nr1)):
    
    fpath = makeFilePath(nr1[i], 'lmp', 'ss', ['thermo'])
    [t, T, pe, msd] = readData(fpath)

    t = [j*dt1[i] for j in t[1:len(t)-1]]
    msd = [j for j in msd[1:len(msd)-1]]
    
    ax1.loglog(t, msd, linestyle = '-', lw = 0.9, label = str(label1[i]))
    
ax1.legend()

nr2 = [6,7]
label2 = [1e-4, 2.5e-3,0.01,0.25,1]
dt2 = [5e-7, 2.5e-6, 5e-6, 2.5e-5, 5e-5]











"""
dt2 = 1e-4
nr2 = [1,2,3,4,5]
label2 = [0.01, 0.05, 0.1, 0.5, 1]

fig2, ax2 = mpl.pyplot.subplots()

for i in range(len(nr2)):
    fname = os.getcwd().removesuffix('/ld_analysis/log') + '/LD/lmp/Data{nr}/thermo.dat'.format(nr = nr2[i])
    fobj = open(fname, mode = "r", encoding = "utf-8")
    fobj.readline()
    
    
    ax2.loglog(t2, msd2, lw = 0.8, label = "$D^*=~$" + str(label2[i]))

ax2.set(xlim=(-10, 1000), ylim = (1e-2,1010), xlabel=r'$t_2^*~\big(\sigma\sqrt{m/k_BT}\big)$', ylabel='$\log{(MSD)}$')
legend = ax2.legend(fontsize=8, title='LAMMPS$~(t_2^*)$')
mpl.pyplot.setp(legend.get_title(), fontsize=8)


#fig2.savefig("msd_lmp.png", dpi=600,bbox_inches='tight')

dt1 = 1e-3
nr1 = 1
t1 = []
msd1 = []
fac=316.23

fname = os.getcwd().removesuffix("/ld_analysis/lmp") + "/LD/cpp/Data{nr1}/relaxa_eq/msd{nr1}.dat".format(nr1=nr1)
f_obj = open(fname, "r", encoding = 'utf-8')
for line in f_obj:
    line = line.split(sep=" ")
    t1.append(liteval(line[0])*fac*dt1)
    msd1.append(liteval(line[1]))


ax1 = ax2.twiny()
ax1.semilogy(t1[0:len(t1)-1:100],msd1[0:len(t1)-1:100],'ko',ms=1,label=r'$m^*=10^{-5}$')
ax1.set(xlim=(-0.1,10), xlabel = r'$t_1^*~(\sigma^2 / D)$')
legend = ax1.legend(fontsize=8,loc='best',bbox_to_anchor=(0.9965,0.39), title='C++$~(t_1^*)$')
mpl.pyplot.setp(legend.get_title(),fontsize=8)


ax1 = ax2.twiny()
ax1.semilogy(t1[0:len(t1)-1:100],msd1[0:len(t1)-1:100],'ko',ms=1,label=r'$m^*=10^{-5}$')
ax1.set(xlim=(-10,1000))
ax1.set_xlabel(r'$t_1^* / \sqrt{m^*} ~\big(\sigma\sqrt{m/k_BT}\big)$', labelpad=12)
legend = ax1.legend(fontsize=8,loc='best',bbox_to_anchor=(0.9965,0.39), title='C++$~(t_1^*)$')
mpl.pyplot.setp(legend.get_title(),fontsize=8)

#fig2.savefig("msd_t2.png", dpi=600,bbox_inches='tight')
"""