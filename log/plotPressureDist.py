#print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "D", "*"]
colors = ["C0", "C1", "C2"]
lw = 0.8
ms = 3
fs = "full"
fontsize = 14

option = "tavg"
if(option == "tavg"):
    
    nr_list = [33]
    lamb_list = [1e-2]
    Pe_list = [1]
    Lx = 200
    Ly = 100
    
    fig1, [ax1, ax2, ax3] = mpl.pyplot.subplots(3,1)
    
    ax1.axvline(x = 100, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw, alpha = 0.6)
    ax2.axvline(x = 100, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw, alpha = 0.6)
    
    for i in range(len(nr_list)):
        
        fpath = "../../LD/LD-cpp/Data{nr}/pressure_tavg_4.dat".format(nr = nr_list[i])
        figpath = "../../LD/LD-cpp/imgs/pressure_tavg.png"
        [x, Pin_xx, Pin_xy, Pin_yx, Pin_yy] = file_utils.readData(fpath, 1)
        
        x = [j-Lx/(2*len(x)) for j in x]
        Pdiff = [(Pin_xx[i] - Pin_yy[i]) for i in range(len(x))]
        
        ax1.plot(x, Pin_xx, marker = markers[0], ms = ms, color = colors[i], ls = "--", 
                lw = lw, fillstyle = fs)
        """
        ax1.plot(x, Pin_xy, marker = markers[1], ms = ms, color = colors[i], ls = "--", 
                lw = lw, fillstyle = fs)
        ax1.plot(x, Pin_yx, marker = markers[2], ms = ms, color = colors[i], ls = "--", 
                lw = lw, fillstyle = fs)
        """
        ax1.plot(x, Pin_yy, marker = markers[3], ms = ms, color = colors[i], ls = "--", 
                lw = lw, fillstyle = fs)
        ax2.plot(x, Pdiff, marker = markers[0], ms = ms, lw = lw)
        #line.set_dashes([1, 1, 10, 2])
        
        print(min(Pdiff), max(Pdiff))
    
    ax1.set(xlim = (0,Lx), aspect = "auto")
    ax1.set_xlabel(r"$x~(\sigma)$", fontsize = fontsize)
    ax1.set_ylabel(r"$\mathcal{P}_{xx}$", fontsize = fontsize)
    
    ax2.set(xlim = (0, Lx), aspect = "auto")
    ax2.set_xlabel(r"$x~(\sigma)$", fontsize = fontsize)
    ax2.set_ylabel(r"$\mathcal{P}_{N} - \mathcal{P}_{T}$", fontsize = fontsize)
    
    fig1.set_figheight(16)
    fig1.set_figwidth(8)
    
    #fig1.savefig(figpath, dpi = 600, bbox_inches = "tight")