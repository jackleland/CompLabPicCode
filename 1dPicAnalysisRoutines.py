# -*- coding: utf-8 -*-

# can't import normally becuase name starts with 1
PIC = __import__('1dPIC')
import IO
import numpy as np
import matplotlib.pyplot as plt

####################################################################
#                        Main Running Script                       #
#################################################################### 


def runNAnalysis(nCell, nPart, N=5, length=4.*np.pi):
    '''
    Run simulation N times with a given nCell & nPart and then return
    the values and standard errors for Omega, Gamma and Noise Amplitude
    '''
    vals = [[],[],[]]
#    plt.figure()
    for i in range(N):
        [tArr, fh, time, a_n, da_n, pt, params, cov, chi] = PIC.setup(2,nPart,nCell,length)
    #    plt.plot(np.pi/pt[1:])
        w2 = np.pi/pt[1]
        a = pt[1:]-pt[0:-1]
        b = a[np.where(a > np.mean(a)*0.8)]
        w = (np.pi)/(np.mean(b))
        dw = np.std(b)/np.sqrt(len(b))
        g = params[1]
        dg = np.sqrt(cov[1][1])
#        print("w = ",w," +/- ",dw)    
#        print("g = ",g," +/- ",dg)
#        print("A_n = ",A_n," +/- ",dA_n)
        vals[0].append(w)
        vals[1].append(g)
        vals[2].append(a_n)
#        data.append(fh)
#        plt.plot(tArr,fh)
    
#    plt.plot(tArr, np.mean(data,0), color='k', lw=1.5)
#    plt.xlabel("Time [Normalised]")
#    plt.ylabel("First harmonic amplitude [Normalised]")
#    plt.yscale('log')
    
    W = np.mean(vals[0])
    dW = standardError(vals[0])
    G = np.mean(vals[1])
    dG = standardError(vals[1])
    A_n = np.mean(vals[2])
    dA_n = standardError(vals[2])
#    print("w = ",W," +/- ",dW)    
#    print("g = ",G," +/- ",dG)
#    print("A_n = ",A_n," +/- ",dA_n)
    return [W, dW, G, dG, A_n, dA_n]



####################################################################
#                       Simulation Timing                          #
####################################################################    

def repeatedNpartRuns(n=20,Ncell=20,scalefactor=5000,plotFl=False):
    data = [[],[]]
    for i in range(n):
        [tArr, fh, time, a_n, da_n, pt, params, cov, chi] = PIC.setup(2,scalefactor*(i+1),Ncell)
        data[0].append(scalefactor*(i+1))
        data[1].append(time)
        
    if plotFl:
        plt.figure()
        plt.plot(data[0],data[1])
        plt.xlabel("$N_{part}$")
        plt.ylabel("Simulation Time (seconds)")  
        plt.title("Simulation Time as a Function of $N_{part}$")
        plt.ioff() # This so that the windows stay open
        plt.show()
    IO.writeCSVFile('dataT-Np-Py2.csv',data)

    
def repeatedNcellRuns(n=20,Npart=10000,scalefactor=5,plotFl=False):
    data = [[],[]]
    for i in range(n):
        [tArr, fh, time, a_n, da_n, pt, params, cov, chi] = PIC.setup(2,Npart,scalefactor*(i+1)) 
        data[0].append(scalefactor*(i+1))
        data[1].append(time)
        
    if plotFl:
        plt.figure()
        plt.plot(data[0],data[1])
        plt.xlabel("$N_{cell}$")
        plt.ylabel("Simulation Time (seconds)")   
        plt.title("Simulation Time as a Function of $N_{cell}$")
        plt.ioff() # This so that the windows stay open
        plt.show()
    IO.writeCSVFile('dataT-Nc-Cy.csv',data)
    
    
    
####################################################################
#                        Parameter Sweeps                          #
####################################################################
   
def runParamaterSweep():
    dataNC=[[],[],[],[],[],[]]
    dataNP=[[],[],[],[],[],[]]
    n_cell = [20,40,60,80,100,120,140,160,180,200,220,240]
    n_part = [1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000]
    for nc in n_cell:
        print("\n********************************************")
        print("N_part = ",10000,", N_cell = ",nc)
        print("********************************************\n")
        [W, dW, G, dG, A_n, dA_n] = runNAnalysis(nc,10000,10)
        dataNC[0].append(W)
        dataNC[1].append(dW)
        dataNC[2].append(G)
        dataNC[3].append(dG)
        dataNC[4].append(A_n)
        dataNC[5].append(dA_n)
#    plt.figure()
#    plt.plot(n_cell,dataNC[0])
        
    # save data to an output csv file
    fileData = ["N_cell",n_cell,"W", dataNC[0], "dW", dataNC[1], "G", dataNC[2], "dG", dataNC[3], "A_n", dataNC[4], "dA_n", dataNC[5]]
    IO.writeCSVFile('dataNC1.csv',fileData)
    
    for nparts in n_part:
        print("\n********************************************")
        print("N_part = ",nparts,", N_cell = ",20)
        print("********************************************\n")

        [W, dW, G, dG, A_n, dA_n] = runNAnalysis(20,nparts,10)
        dataNP[0].append(W)
        dataNP[1].append(dW)
        dataNP[2].append(G)
        dataNP[3].append(dG)
        dataNP[4].append(A_n)
        dataNP[5].append(dA_n)
#    plt.figure()
#    plt.plot(n_cell,dataNC[0])
        
    # save data to an output csv file
    fileData = ["N_part",n_part,"W", dataNP[0], "dW", dataNP[1], "G", dataNP[2], "dG", dataNP[3], "A_n", dataNP[4], "dA_n", dataNP[5]]
    IO.writeCSVFile('dataNP1.csv',fileData)
        
    
def runLengthSweep():
    dataL = [[],[],[],[],[],[]]
    lVals = np.linspace(4,10,10)*np.pi
    for l in lVals:
        print("\n********************************************")
        print("N_part = ",100000,", N_cell = ",20,", l = ",l)
        print("********************************************\n")
        [W, dW, G, dG, A_n, dA_n] = runNAnalysis(20,200000,5,l)
        dataL[0].append(W)
        dataL[1].append(dW)
        dataL[2].append(G)
        dataL[3].append(dG)
        dataL[4].append(A_n)
        dataL[5].append(dA_n)
    
    # save data to an output csv file
    IO.writeCSVFile('dataL1.csv',dataL)
    
    
if __name__ == "__main__":
    repeatedNpartRuns()
