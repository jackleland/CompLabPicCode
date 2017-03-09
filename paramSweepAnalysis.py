# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import IO
from scipy.optimize import curve_fit 
import ffunctions as ff


nc_data = IO.readCSVFile('dataNC1.csv')
np_data = IO.readCSVFile('dataNP1.csv')

#####################################################################
##                        Number of Cells                           #
#####################################################################
#
## omega measurements
#plt.subplot(2,3,1)
#plt.errorbar(nc_data[0],nc_data[1],nc_data[2])
#plt.xlabel(r'Number of cells ($N_{cell}$)')
#plt.ylabel(r'Mean frequency of first harmonic ($\bar{\omega}$)')
#
## gamma measurements
#plt.subplot(2,3,2)
#plt.errorbar(nc_data[0],nc_data[3],nc_data[4])
#plt.xlabel(r'Number of cells ($N_{cell}$)')
#plt.ylabel(r'Mean landau damping rate ($\bar{\gamma}$)')
#
## Noise measurements
#plt.subplot(2,3,3)
#plt.errorbar(nc_data[0],nc_data[5],nc_data[6])
#plt.xlabel(r'Number of cells ($N_{cell}$)')
#plt.ylabel(r'Mean noise amplitude ($\bar{A}_N$)')
#
#####################################################################
##                      Number of Particles                         #
#####################################################################
#
## omega measurements
#plt.subplot(2,3,4)
#plt.errorbar(np_data[0],np_data[1],np_data[2],color='g')
#plt.xlabel(r'Number of particles ($N_{part}$)')
#plt.ylabel(r'Mean frequency of first harmonic ($\bar{\omega}$)')
#plt.xscale('log')
#plt.xlim(8.5e2,1.5e6)
#
## gamma measurements
#plt.subplot(2,3,5)
#plt.errorbar(np_data[0],np_data[3],np_data[4],color='g')
#plt.xlabel(r'Number of particles ($N_{part}$)')
#plt.ylabel(r'Mean landau damping rate ($\bar{\gamma}$)')
#plt.xscale('log')
#plt.xlim(8.5e2,1.5e6)
#
## Noise measurements
#plt.subplot(2,3,6)
#plt.errorbar(np_data[0],np_data[5],np_data[6],color='g')
#plt.xlabel(r'Number of particles ($N_{part}$)')
#plt.ylabel(r'Mean noise amplitude ($\bar{A}_N$)')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(8.5e2,1.5e6)

N = np_data[0]
A = np_data[5]
dA = np_data[6]

snr = 1.0/A
dsnr = dA/(A)**2

logN = np.log(np_data[0])
logA = np.log(np_data[5])
dlogA = np.divide(np_data[6],np_data[5])

guess = [-2.65,15]
popt,pcov = curve_fit(ff.linear, logN[0:5], logA[0:5], p0=guess, sigma=dlogA[0:5])
print(popt)
fit = ff.linear(logN[0:5],*popt)

guess1 = [1.0,1.0]
popt1,pcov1 = curve_fit(ff.squareRoot, N[0:6], snr[0:6], p0=guess1, sigma=dsnr[0:6])
print(popt1)
fit1 = ff.squareRoot(N[0:6],*popt)

plt.errorbar(logN,logA,dlogA)
#plt.plot(np_data[0][0:6],np_data[5][0:6])
plt.plot(logN[0:6],fit)
plt.xlabel(r'Log(Number of particles) ($N_{part}$)')
plt.ylabel(r'Log(Mean noise amplitude) ($\bar{A}_N$)')
#plt.xscale('log')

plt.figure()
plt.errorbar(N[0:6],snr[0:6],dsnr[0:6])
plt.plot(N[0:6],fit1)

