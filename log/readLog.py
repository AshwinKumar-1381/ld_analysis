"""
Extracts columns of thermo data from log.lammps and writes to thermo.dat

Author          : Ashwin Kumar M
Date created    : 16.10.24
Last modified   : 16.10.24

-------- INPUTS -------
argv[1] = nr    (file no.)
argv[2] = nRuns (no. of run sets present in log.lammps)
-----------------------
"""

print("\033[J\033[H",end='') # Clear screen

import os
import sys

def line2array(line):
    line = line.removesuffix("\n").split(sep = " ")
    line2 = []
    for element in line:
        if(element != ''): line2.append(element)
    if(line2 == []): line2.append('0')
    return(line2)

def array2line(array):
    line = ""
    for i in range(len(array)):
        line += str(array[i]) + " "
    return(line.removesuffix(" "))

args = sys.argv

fname1 = os.getcwd().removesuffix('/ld_analysis/log') + \
        '/LD/lmp/Data{nr}/log.lammps'.format(nr = int(args[1]))
fobj1 = open(fname1, mode = "r", encoding = "utf-8")

fname2 = fname1.removesuffix("/log.lammps") + "/thermo.dat"
fobj2 = open(fname2, mode = "w")

run_nr = 0
for line in fobj1:
    array = line2array(line)
    
    if(run_nr == int(args[2])):
        if(array[0] == 'Loop'): break
        else:
            line = array2line(array) 
            fobj2.write(line + "\n")
    
    if(array[0] == 'Per'): 
        run_nr += 1

fobj1.close()
fobj2.close()
