"""
 Script to export columns of data from dump/thermo.dat file to plot potential energies
 Author         : Ashwin Kumar M (CH23S006)
 Date created   : 13.08.24
 Last modified  : 16.10.24
"""

print("\033[J\033[H",end='') # Clear screen

from ast import literal_eval as liteval
import matplotlib as mpl
import os

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

nr = [11]
dt = 5e-5
type1 = 'lmp'
type2 = 'ss'
plotInterval = 10

fig, ax = mpl.pyplot.subplots()

for j in range(len(nr)):

    logPath = makeFilePath(nr[j], type1, type2, ['log'])
    
    if(type1 == 'cpp'):
        PeA = readLogParam('Pe_A', logPath)
        PeB = readLogParam('Pe_B', logPath)
        dt = readLogParam('dimensionless_mass', logPath)
        N = readLogParam('N', logPath)
        
        dumpPath = makeFilePath(nr[j], type1, type2, ['dump'])
        
        if(type2 == 'eq'): 
            [t, pe] = readData(dumpPath)
            lent = len(t)
            for i in range(lent): t[i] *= dt
            
            ax.plot(t[0 : lent - 1 : plotInterval], pe[0 : lent - 1 : plotInterval],
                    lw = 0.8, label = "total PE")
        
        elif(type2 == 'ss'): 
            [t, pe, peA, peB] = readData(dumpPath)
            lent = len(t)
            for i in range(lent): t[i] *= dt
            
            ax.plot(t[0 : lent - 1 : plotInterval], pe[0 : lent - 1 : plotInterval],
                    lw = 0.8, label = "total PE")
            ax.plot(t[0 : lent - 1 : plotInterval], peA[0 : lent - 1 : plotInterval],
                    lw = 0.8, label = "partivle A PE")
            ax.plot(t[0 : lent - 1 : plotInterval], peB[0 : lent - 1 : plotInterval],
                    lw = 0.8, label = "particle B PE")
        
    elif(type1 == 'lmp'):
        N = readLogParam('N', logPath)
        thermoPath = makeFilePath(nr[j], type1, type2, ['thermo'])
        [t, T, pe] = readData(thermoPath)
        
        lent = len(t)
        for i in range(lent): 
            t[i] *= dt
            pe[i] *= N
        
        ax.plot(t[0 : lent - 1 : plotInterval], pe[0 : lent - 1 : plotInterval],
                lw = 0.8, label = "total PE")

        figPath = makeFilePath(nr[j], type1, type2, ['fig'])
    
ax.set_xlabel("$t^*$", fontsize = 15)
ax.set_ylabel("Potential Energy", fontsize = 12)
ax.legend()

#fig.savefig(figPath, dpi = 600, bbox_inches = 'tight')