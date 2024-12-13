print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "*", "D"]
lw = 0.5
ms = 3
fontsize = 10
fac = 1

option = "tavg"
if(option == "frames"):
    # Params for frames-based data
    nr = 21
    Lx = 200
    Pe = "1.0"
    lamb = " 10^{-3}" 
    pdfSamples = [11e3, 12e3, 13e3, 14e3, 15e3]
    nSamples = len(pdfSamples)
    sampleNum = [0, 1, 2, 3, 4]
    
    fpath = "../../LD/LD-cpp/Data{nr}/pdf_frames.dat".format(nr = nr)
    figpath = "../../LD/LD-cpp/imgs/pdf_frames_{nr}.png".format(nr = nr)
    pdfData = file_utils.readData(fpath, nSamples)
    
    fig, [axA, axB] = mpl.pyplot.subplots(1, 2)
    for i in sampleNum:
        x = pdfData[3*i + 0]
        pdfA = pdfData[3*i + 1]
        pdfB = pdfData[3*i + 2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        
        axA.plot(x, pdfA, "b-", marker = markers[i], lw = lw, ms = ms, label = str(pdfSamples[i]))
        axB.plot(x, pdfB, "r-", marker = markers[i], lw = lw, ms = ms, label = str(pdfSamples[i]))
    
    axA.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axB.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axA.set_xlabel(r"$x$", fontsize = fontsize)
    axB.set_xlabel(r"$x$", fontsize = fontsize)
    axA.set_ylabel(r"$\phi_{A} (x), ~\phi_{B} (x)$", fontsize = fontsize)
    axA.legend(title = "time", fontsize = fontsize - 3)
    axB.legend(title = "time", fontsize = fontsize - 3)
    fig.suptitle(r"Simulation Parameters : $Pe_s = {Pe},~ \lambda = {lamb}$".format(Pe=Pe,lamb=lamb), 
                 fontsize = fontsize)
    
    #fig.savefig(figpath, dpi = 600, bbox_inches = "tight")
    
elif(option == "tavg"):
    # Params for time-averaged data
    nr_list = [21]
    simParams = [1,5]
    lamb = "10^{-3}"
    Lx = 200.0
    
    fig, [axA, axB] = mpl.pyplot.subplots(1, 2)
    for i in range(len(nr_list)):
        fpath = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = nr_list[i])
        figpath = "../../LD/LD-cpp/imgs/pdf_tavg_lamb_1e-2.png"
        
        pdfData = file_utils.readData(fpath, 1)
        
        x = pdfData[0]
        pdfA = pdfData[1]
        pdfB = pdfData[2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        
        axA.plot(x, pdfA, "b-", marker = markers[i], lw = lw, ms = ms, label = str(simParams[i]))
        axB.plot(x, pdfB, "r-", marker = markers[i], lw = lw, ms = ms, label = str(simParams[i]))
    
    axA.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axB.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axA.set_xlabel(r"$x ~(\sigma)$", fontsize = fontsize)
    axB.set_xlabel(r"$x ~(\sigma)$", fontsize = fontsize)
    axA.set_ylabel(r"$<\phi_{A}>(x), ~<\phi_{B}>(x)$", fontsize = fontsize)
    axA.legend(title = r"$Pe_s$", fontsize = fontsize - 4)
    axB.legend(title = r"$Pe_s$", fontsize = fontsize - 4)
    fig.suptitle(r"Simulation Parameters : $\lambda = {lamb}$".format(lamb=lamb), fontsize = fontsize)
    
    #fig.savefig(figpath, dpi = 600, bbox_inches = "tight")    