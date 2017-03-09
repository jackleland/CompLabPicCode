# CompLabPicCode
PIC code from first year CDT computational labs

Main Files:
 - 1dPIC.py - main simulation script
 - 1dPicAnalysisRoutines.py - functions for running 1dPIC in a variety of configurations

Analysis/Plotting Files:
 - paramSweepAnalysis.py
 - BoxLengthAnalysis.py
 - tSimAnalysis.py
 
Cython Files:
 - cOptimisation.pyx - original algorithm in cython
 - setup.py - cython make file
 - cOptimisation.c - c object of algorithm, callable in python

Additional Modules:
 - IO.py - module for reading and writing csv files using numpy arrays
 - ffunctions.py - assorted functions for fitting with scipy's curve_fit() function
