# -*- coding: utf-8 -*-
import numpy as np

# This function will generate a perfect Gaussian 
def gaussian(x, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    y0 = params[3]
    
    return y0 + A*np.exp(-(x-x0)**2/(2*c*c))
    
def gaussian(x, y0, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    
    return y0 + A*np.exp(-(x-x0)**2/(2*c*c))
    
# This function will generate a perfect super gaussian of power n    
def superGaussian(x, n, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    y0 = params[3]

    return y0 + A*np.exp(-((x-x0)/(2*c))**n)
    
# This function will generate a bi-gaussian with one being a power 6 super gaussian    
def biGaussian(x, *params):
    params1 = [params[0],params[1],params[2],params[3]]
    params2 = [params[4],params[5],params[6],params[3]]
    
    return gaussian(x, *params1) + superGaussian6(x, *params2)
    
def polyGaussian(x, n, *params):
    y0 = params[0]
    output = np.zeros(np.size(x))
    for i in range(0,n):
        A = params[(3*i)+1]
        x0 = params[(3*i)+2]
        c = params[(3*i)+3]
        output = output + A*np.exp(-((x-x0)/(2*c))**2)
        
    return y0 + output
        
    
# This function will generate a bi-gaussian with one being a power 6 super gaussian
# using an alternate method for calculation    
def biGaussianAlt(x, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    y0 = params[3]
    A2 = params[4]
    x02 = params[5]
    c2 = params[6]
    
    return y0 + A*np.exp(-((x-x0)/(2*c))**2) + A2*np.exp(-((x-x02)/(2*c2))**6)
    
     
# This function will generate a power 4 super gaussian by calling the superGaussian function    
def superGaussian4(x, *params):
    return superGaussian(x,4,*params)
    
# This function will generate a power 6 super gaussian by calling the superGaussian function    
def superGaussian6(x, *params):
    return superGaussian(x,6,*params)
    
# This function will generate a power 8 super gaussian by calling the superGaussian function    
def superGaussian8(x, *params):
    return superGaussian(x,8,*params)
    
def polyGaussian3(x, *params):
    return polyGaussian(x,3,*params)

def polyGaussian4(x, *params):
    return polyGaussian(x,4,*params)

def polyGaussian5(x, *params):
    return polyGaussian(x,5,*params)
    
def polyGaussian6(x, *params):
    return polyGaussian(x,6,*params)
    
def exponential(x, *params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*np.exp(b*x) + c
    
def linear(x, *params):
    m = params[0]
    c = params[1]
    return m*x + c
    
def quadratic(x, *params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*(x**2) + b*x + c
    
def cubic(x, *params):
    a = params[0]
    b = params[1]
    c = params[2]
    d = params[3]    
    return a*(x**3) + b*(x**2) + c*x + d
    
def cubnoquad(x, *params):
    a = params[0]
    b = params[1]
    c = params[2]    
    return a*(x**3) + b*x + c 
    
def sinusoid(x, *params):
    a = params[0]
    b = params[1]
    c = params[2]
    d = params[3]
    return a*np.sin(b*(x + c)) + d
