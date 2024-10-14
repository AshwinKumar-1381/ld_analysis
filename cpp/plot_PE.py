"""
 Script to export columns of data from dump file to .dat format
 Author: Ashwin Kumar M (CH23S006)
 Date created: 13.08.24
 Last modified: 24.09.24
"""

print("\033[J\033[H",end='') # Clear screen

from ast import literal_eval as liteval
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import matplotlib as mpl
import os

mpl.pyplot.rc('text', usetex=True)

nr = [1]
restart_nr = 0
type = 'eq'
plot_interval = 10

for j in range(len(nr)):
    # log file path
    fpath = os.getcwd().removesuffix("code") + "LD/Data{nr}/log.dat".format(nr = nr[j])
    f_obj = open(fpath, mode = "r", encoding = "utf-8")
    for line in f_obj:
        line = line.removesuffix("\n").split(sep = " ")
        if(line[0] == "Pe_A"): PeA = liteval(line[1])
        if(line[0] == "Pe_B"): PeB = liteval(line[1])
        if(line[0] == "dimensionless_mass"): dt = liteval(line[1])
        if(line[0] == "N"): N = liteval(line[1])

    f_obj.close()

    # dump.dat file path
    fpath = os.getcwd().removesuffix('code') + 'LD/Data' + str(nr[j]) + '/relaxa_' +\
            type + '/dump_' + type + '.dat'

    i = 0
    time = []
    pe_total = []
    pe_A = []
    pe_B = []

    # energy plot for equilibrium run
    if(type == 'eq'):
    
        # open dump.dat file
        f_obj = open(fpath, mode = 'r', encoding = 'utf-8')
        for line in f_obj: 
            line = line.split(sep = ' ')
            if(i==0): print('')
            elif(i%plot_interval == 0):
                time.append(liteval(line[0])*(dt/10))
                pe_total.append(liteval(line[1]))
            i += 1

        fig, ax = mpl.pyplot.subplots()

        ax.plot(time, pe_total, '-', lw = '0.5 ')

    # energy plot for steady-state run
    else: 
        if(restart_nr != 0): fpath = fpath.replace("relaxa_ss", "restart_" + str(restart_nr))
    
        # open dump.dat file
        f_obj = open(fpath, mode = 'r', encoding = 'utf-8')
    
        for line in f_obj: 
            line = line.split(sep = ' ')
            if(i==0): print('')
            elif(i%plot_interval == 0):
                time.append(liteval(line[0])*dt)
                pe_total.append(liteval(line[1]))
                pe_A.append(liteval(line[2]))
                pe_B.append(liteval(line[3]))
            i += 1
        
        fig, ax = mpl.pyplot.subplots()
        ax.plot(time, pe_total, '-', ms=0.5, lw = '0.5 ', color = 'green', label=r'E$_A$ + E$_B$')
        ax.plot(time, pe_A, '-', ms=0.5, lw = '0.5', color = '#2D51F5', label = r'E$_A$')
        ax.plot(time, pe_B, '-', ms=0.5, lw = '0.5', color = '#F34943', label = r'E$_B$')
        ax.legend()
        ax.set_title("Hard sphere potential energy " + r'$(N = {N}, Pe_A = {PeA}, ~Pe_B = {PeB})$'\
                     .format(N = N, PeA = PeA, PeB = PeB), fontsize = 10)

ax.set(xlabel = r'time ($\tau$)', ylabel = r'E$_{pot} ~ (k_B T)$')

#fpath = fpath.removesuffix('/dump_' + type + '.dat') + '/pe_ss_{NR}'.format(NR = nr[j]) + '.png'
#fig.savefig(fpath, dpi = 500, bbox_inches = 'tight')