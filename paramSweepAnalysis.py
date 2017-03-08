# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import csv

def readCSVFile(filename):
    f = open(filename, 'r')
    reader = csv.reader(f)
    data = []
    
    while True:
        try:
            d = reader.next()[0]
        except:
            break
        try:
            data.append(np.array(d.strip("[]").split(", ")).astype(np.float))
        except:
            continue
    
    f.close()
    return data

nc_data = readCSVFile('dataNC1.csv')
np_data = readCSVFile('dataNP1.csv')

####################################################################
#                        Number of Cells                           #
####################################################################

# omega measurements
fig0 = plt.figure(0)
plt.errorbar(nc_data[0],nc_data[1],nc_data[2])
plt.xlabel(r'Number of cells ($N_{cell}$)')
plt.ylabel(r'Mean frequency of first harmonic ($\bar{\omega}$)')

# gamma measurements
plt.figure(1)
plt.errorbar(nc_data[0],nc_data[3],nc_data[4])
plt.xlabel(r'Number of cells ($N_{cell}$)')
plt.ylabel(r'Mean landau damping rate ($\bar{\gamma}$)')

# Noise measurements
plt.figure(2)
plt.errorbar(nc_data[0],nc_data[5],nc_data[6])
plt.xlabel(r'Number of cells ($N_{cell}$)')
plt.ylabel(r'Mean noise amplitude ($\bar{A}_N$)')

####################################################################
#                      Number of Particles                         #
####################################################################

# omega measurements
plt.figure(3)
plt.errorbar(np_data[0],np_data[1],np_data[2],color='g')
plt.xlabel(r'Number of particles ($N_{part}$)')
plt.ylabel(r'Mean frequency of first harmonic ($\bar{\omega}$)')
plt.xscale('log')
plt.xlim(8.5e2,1.5e6)

# gamma measurements
plt.figure(4)
plt.errorbar(np_data[0],np_data[3],np_data[4],color='g')
plt.xlabel(r'Number of particles ($N_{part}$)')
plt.ylabel(r'Mean landau damping rate ($\bar{\gamma}$)')
plt.xscale('log')
plt.xlim(8.5e2,1.5e6)


# Noise measurements
plt.figure(5)
plt.errorbar(np_data[0],np_data[5],np_data[6],color='g')
plt.xlabel(r'Number of particles ($N_{part}$)')
plt.ylabel(r'Mean noise amplitude ($\bar{A}_N$)')
plt.xscale('log')
plt.xlim(8.5e2,1.5e6)

