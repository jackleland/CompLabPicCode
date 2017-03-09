# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import IO
from scipy.optimize import curve_fit 
import ffunctions as ff


lVals = np.linspace(2,10,5)/2
L_data = IO.readCSVFile('dataL.csv')



####################################################################
#                           Box Length                             #
####################################################################

# omega measurements
plt.subplot(1,2,1)
plt.errorbar(lVals,L_data[0],L_data[1])
plt.xlabel(r'Length of box ($\hat{L}$)')
plt.ylabel(r'Mean frequency of first harmonic ($\bar{\omega}$)')

# gamma measurements
plt.subplot(1,2,2)
plt.errorbar(lVals,L_data[2],L_data[3])
plt.xlabel(r'Length of box ($\hat{L}$)')
plt.ylabel(r'Mean landau damping rate ($\bar{\gamma}$)')

## Noise measurements
#plt.subplot(1,3,3)
#plt.errorbar(lVals,L_data[4],L_data[5])
#plt.xlabel(r'Number of cells ($N_{cell}$)')
#plt.ylabel(r'Mean noise amplitude ($\bar{A}_N$)')


