# -*- coding: utf-8 -*-
'''
Module for reading and writing numpy arrays to csv files easily
'''

import csv
import numpy as np

def readCSVFile(fileName):
    '''
    Reads rows of data from a CSV (fileName) into a numpy array and returns it
    Note: will skip any rows that cannot be converted from the saved numpy 
    array format e.g. strings/labels.
    '''
    f = open(fileName, 'r')
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
    
    
def writeCSVFile(fileName,fileData):
    '''
    Writes numpy array data (fileData) row by row to a csv with name fileName
    '''
    f = open(fileName, 'w') 
    writer = csv.writer(f)
    writer.writerows([[x] for x in fileData])
    f.close()

