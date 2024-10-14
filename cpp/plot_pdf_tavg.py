"""
 Code to read time-averaged Particle Density Function (PDF)
 data corresponding to single and multiple runs
 Author: Ashwin Kumar M (CH23S006)
 Date created: 09.09.24
 Last modified: 26.09.24
"""

print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
from ast import literal_eval as liteval
import os
import numpy as np

mpl.pyplot.rc('text', usetex=True)

nr_list = [16]
markers = ["o", '^', "s", "*"]
line_style = ["-", "-", "-", "-", "--"]

lw = 0.5
cblue = "#2D51F5"
cred = "#F34943"

if(len(nr_list) == 1):
    fig, ax = mpl.pyplot.subplots(figsize = (5,5))
    x = []
    pdfA = []
    pdfB = []
    fpath = os.getcwd().removesuffix("/code") + \
            "/LD/Data{nr}/relaxa_ss/pdf_{nr}.dat".format(nr = nr_list[0])
    
    f_obj = open(fpath, mode = "r", encoding = "utf-8")
    
    line = f_obj.readline().split(sep = " ")
    n_bin = liteval(line[0])
    L = liteval(line[1])
    N = liteval(line[2])
    PeA = liteval(line[3])
    PeB = liteval(line[4].removesuffix("\n"))
    
    f_obj.readline().split(sep = " ")

    for j in range(n_bin):
        line = f_obj.readline().split(sep = " ")
        x.append(liteval(line[0])*L/n_bin)
        pdfA.append(liteval(line[1]))
        pdfB.append(liteval(line[2].removesuffix("\n")))
    
    # Plot PDf of particles A and B vs x
    ax.plot(x, pdfA, ls = line_style[0], lw = lw, color = cblue, marker = markers[0], ms = 3,
            mec = cblue, mfc = "white", label = r"$Pe_A = $" + str(PeA))
    ax.plot(x, pdfB, ls = line_style[0], lw = lw, color = cred, marker = markers[0], ms = 3,
            mec = cred, mfc = "white", label = r"$Pe_B = $" + str(PeB))
    ax.legend(fontsize = 8)    
    ax.set(xlabel = r"$x ~(\sigma)$", ylabel = r"$\rho (x)$")
    
    # Compute the integral of the PDf
    sum = [0]
    for i in range(len(x) - 1):
        sum.append(sum[i] + 0.5*(x[i+1] - x[i])*(pdfB[i+1] + pdfB[i]))
    for i in range(len(sum)):
        sum[i] *= 2*L/N
        sum[i] = round(sum[i],5)
    
    # Plot the integral of the PDF vs x
    fig2, ax2 = mpl.pyplot.subplots(figsize = (5,5))
    ax2.axhline(y = 0, xmin = 0, xmax = 1, ls = "--", c = 'k', lw = 0.8, alpha = 0.7)
    ax2.axhline(y = 1, xmin = 0, xmax = 1, ls = "--", c = 'k', lw = 0.8, alpha = 0.7)
    ax2.axvline(x = 40, ymin = 0, ymax = 1, ls = "--", c = 'k', lw = 0.8, alpha = 0.7)
    ax2.axvline(x = 65, ymin = 0, ymax = 1, ls = "--", c = 'k', lw = 0.8, alpha = 0.7)
    ax2.plot(x, sum, ls = line_style[0], lw = lw, color = cred, marker = markers[0], ms = 3,
             mec = cred, mfc = 'white', label = r"$Pe_B = $" + str(PeB))
    ax2.set(xlabel = r"$x ~(\sigma)$", ylabel = r"$I (x)$")
    
    fig.savefig(fpath.removesuffix("{nr}.dat".format(nr = nr_list[0])) + "plot_{nr}_{n}.png".format(nr = nr_list[0], n = N),
                dpi = 600, bbox_inches = "tight")
    #fig2.savefig(fpath.removesuffix("{nr}.dat".format(nr = nr_list[0])) + "int_{nr}.png".format(nr = nr_list[0]), dpi = 600, bbox_inches = "tight")
    
else:
    
    fig, [axA, axB] = mpl.pyplot.subplots(1, 2)
    x = []
    pdfA = []
    pdfB = []
    for i in range(len(nr_list)):
        fpath = os.getcwd().removesuffix("/code") + \
                "/LD/Data{nr}/relaxa_ss/pdf_{nr}.dat".format(nr = nr_list[i])
    
        f_obj = open(fpath, mode = "r", encoding = "utf-8")
    
        line = f_obj.readline().split(sep = " ")
        n_bin = liteval(line[0])
        L = liteval(line[1])
        N = liteval(line[2])
        PeA = liteval(line[3])
        PeB = liteval(line[4].removesuffix("\n"))
    
        f_obj.readline().split(sep = " ")
    
        for j in range(n_bin):
            line = f_obj.readline().split(sep = " ")
            x.append(liteval(line[0])*L/n_bin)
            pdfA.append(liteval(line[1]))
            pdfB.append(liteval(line[2].removesuffix("\n")))
    
        axA.plot(x, pdfA, ls = line_style[i], lw = lw, color = cblue, marker = markers[i], ms = 3, 
             mec = cblue, mfc = "white", label = str(PeA))
        axB.plot(x, pdfB, ls = line_style[i], lw = lw, color = cred, marker = markers[i], ms = 3, 
             mec = cred, mfc = "white", label = str(PeB))
    
        axA.set_xlabel(r"x ($\sigma$)")
        axA.set_ylabel(r"$\rho (x)$")
        axA.set_title("Particle A", fontsize = 10)
        axA.legend(title = "$Pe_A$", fontsize = 8)
        axB.set_xlabel("x ($\sigma$)")
        axB.set_title("Particle B", fontsize = 10)
        axB.legend(title = "$Pe_B$", fontsize = 8)
    
        x = []
        pdfA = []
        pdfB = []
    
        f_obj.close()
    #fig.savefig("pdf_plot.png", dpi = 1000, bbox_inches = 'tight')
