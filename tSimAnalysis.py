# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import IO
from scipy.optimize import curve_fit 
import ffunctions as ff


TCy_data = IO.readCSVFile('dataT-Np-Cy.csv')
TPy_data = IO.readCSVFile('dataT-Np-Py.csv')
TPy1_data = IO.readCSVFile('dataT-Np-Py1.csv')
TPy2_data = IO.readCSVFile('dataT-Np-Py2.csv')


#print("cy: ",TCy_data)
#print("py: ",TPy_data)


####################################################################
#                           Box Length                             #
####################################################################

# time measurements
plt.figure
plt.plot(TCy_data[0],TCy_data[1], label='Cython Alg.')
plt.plot(TPy_data[0],TPy_data[1], label='np.where Alg.')
plt.plot(TPy2_data[0],TPy2_data[1], label='np.bincount Alg.')
plt.plot(TPy1_data[0],TPy1_data[1], label='Original Alg.')
plt.xlabel(r'Number of Particles ($N_{part}$)')
plt.ylabel(r'Simulation Time (t)')
plt.legend(loc='best')



