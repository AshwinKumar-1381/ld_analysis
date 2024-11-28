import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m

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
                                      "/imgs/msd_{type1}_{nr}.png".format(type1=type1, nr = extras[1])
    elif(extras[0] == 'msd' and type1 == 'lmp'): fpath += "/msd.dat"
    elif(extras[0] == 'msd' and type1 == 'cpp'): fpath += "/relaxa_{type2}/msd.dat".format(type2 = type2) 
    return(fpath)

def readData(fpath):
    fobj = open(fpath, mode = 'r', encoding = "utf-8")
    nCols = len(fobj.readline().removesuffix("\n").split(sep = " "))
    Cols = [[] for i in range(nCols)]
    
    for line in fobj:
        line = line.removesuffix("\n").split(sep = " ")
        for i in range(nCols): 
            if(line[i].isalpha() == True): Cols[i].append(line[i])
            else: Cols[i].append(liteval(line[i]))
    
    fobj.close()
    return(Cols)