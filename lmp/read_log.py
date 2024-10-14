print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval

mpl.pyplot.rc('text', usetex=True)

dt2 = 1e-4
nr2 = [2,1,3,4,5]
label2 = [1e-3, 3e-3, 1e-2, 0.1, 1]

fig2, ax2 = mpl.pyplot.subplots()

for i in range(len(nr2)):
    fname = os.getcwd().removesuffix('/ld_analysis/lmp') + '/LD/lmp/Data{nr}/log{nr}.dat'.format(nr = nr2[i])
    f_obj = open(fname, mode = "r", encoding = "utf-8")

    ch = 'n'
    line_nr = -1
    t2 = []
    msd2 = []
    for line in f_obj:
        vals = []
        line = line.removesuffix('\n').split(sep=' ')
    
        if(line[0] == 'Loop'):
            ch = 'n'
            line_nr = -1
    
        if(ch == 'y' and line_nr > 0):
            for element in line:
                if(element != ''): vals.append(liteval(element))
            t2.append(round(vals[0]*dt2,4))
            msd2.append(vals[3])

        if(line[0] == 'Per' or line_nr >= 0): 
            ch = 'y' 
            line_nr += 1
    ax2.semilogy(t2, msd2, lw = 0.8, label = "$D^*=~$" + str(label2[i]))

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

"""
ax1 = ax2.twiny()
ax1.semilogy(t1[0:len(t1)-1:100],msd1[0:len(t1)-1:100],'ko',ms=1,label=r'$m^*=10^{-5}$')
ax1.set(xlim=(-0.1,10), xlabel = r'$t_1^*~(\sigma^2 / D)$')
legend = ax1.legend(fontsize=8,loc='best',bbox_to_anchor=(0.9965,0.39), title='C++$~(t_1^*)$')
mpl.pyplot.setp(legend.get_title(),fontsize=8)
"""

ax1 = ax2.twiny()
ax1.semilogy(t1[0:len(t1)-1:100],msd1[0:len(t1)-1:100],'ko',ms=1,label=r'$m^*=10^{-5}$')
ax1.set(xlim=(-10,1000))
ax1.set_xlabel(r'$t_1^* / \sqrt{m^*} ~\big(\sigma\sqrt{m/k_BT}\big)$', labelpad=12)
legend = ax1.legend(fontsize=8,loc='best',bbox_to_anchor=(0.9965,0.39), title='C++$~(t_1^*)$')
mpl.pyplot.setp(legend.get_title(),fontsize=8)

fig2.savefig("msd_t2.png", dpi=600,bbox_inches='tight')