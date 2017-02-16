#!/usr/bin/env python
#
# Electrostatic PIC code in a 1D cyclic domain

import numpy as np
import time as t
import scipy.signal as spy
import ffunctions as ff

import matplotlib.pyplot as plt # Matplotlib plotting library

try:
    import matplotlib.gridspec as gridspec  # For plot layout grid
    got_gridspec = True
except:
    got_gridspec = False

# Need an FFT routine, either from SciPy or NumPy
try:
    from scipy.fftpack import fft, ifft
except:
    # No SciPy FFT routine. Import NumPy routine instead
    from numpy.fft import fft, ifft

def firstMin(arr):
    mint = arr[0]
    for i in range(len(arr)):
        if mint < arr[i]:
            return i - 1
        else:
            mint = arr[i]
        

def rk4step(f, y0, dt, args=()):
    """ Takes a single step using RK4 method """
    k1 = f(y0, *args)
    k2 = f(y0 + 0.5*dt*k1, *args)
    k3 = f(y0 + 0.5*dt*k2, *args)
    k4 = f(y0 + dt*k3, *args)

    return y0 + (k1 + 2.*k2 + 2.*k3 + k4)*dt / 6.

def calc_density(position, ncells, L):
    """ Calculate charge density given particle positions
    
    Input
      position  - Array of positions, one for each particle
                  assumed to be between 0 and L
      ncells    - Number of cells
      L         - Length of the domain

    Output
      density   - contains 1 if evenly distributed
    """
    # This is a crude method and could be made more efficient
    
    density = np.zeros([ncells])
    nparticles = len(position)
    
    dx = L / ncells       # Uniform cell spacing

    pos = np.divide(position, dx)
    
    for i in range(ncells):
        j = (i+1)%ncells
        sumrange = np.where( (pos >= i) & (pos < i+1) )
        density[i] = density[i] + sum(pos[sumrange] ) - i*len(sumrange[0])
        density[j] = density[j] + (i+1)*len(sumrange[0]) - sum(pos[sumrange] )
#    cellnum = floor(pos).astype(int)
#    relpos = pos - cellnum
#    relposinv = 1 - relpos
#    
#    # counts the number of particles in each cell and allows them to be summed over
#    # as the pos list is sorted. This falls down if either end cell is empty.
#    counts = bincount(cellnum)
#        
#    p_hi = 0
#    p_lo = 0
#    
#    for i in range(ncells):
#        p_lo = p_hi
#        p_hi = p_hi + counts[i]
#        # not using += because it is slower
#        density[i] = density[i] + sum(relposinv[p_lo:p_hi])
#        density[(i+1)%ncells] = density[(i+1)%ncells] + sum(relpos[p_lo:p_hi])
    

        
    # nparticles now distributed amongst ncells
    density *= float(ncells) / float(nparticles)  # Make average density equal to 1
    return density

def periodic_interp(y, x):
    """
    Linear interpolation of a periodic array y at index x
    
    Input

    y - Array of values to be interpolated
    x - Index where result required. Can be an array of values
    
    Output
    
    y[x] with non-integer x
    """
    ny = len(y)
    if len(x) > 1:
        y = np.array(y) # Make sure it's a NumPy array for array indexing
    xl = np.floor(x).astype(int) # Left index
    dx = x - xl
    xl = ((xl % ny) + ny) % ny  # Ensures between 0 and ny-1 inclusive
    return y[xl]*(1. - dx) + y[(xl+1)%ny]*dx

def fft_integrate(y):
    """ Integrate a periodic function using FFTs
    """
    n = len(y) # Get the length of y
    
    f = fft(y) # Take FFT
    # Result is in standard layout with positive frequencies first then negative
    # n even: [ f(0), f(1), ... f(n/2), f(1-n/2) ... f(-1) ]
    # n odd:  [ f(0), f(1), ... f((n-1)/2), f(-(n-1)/2) ... f(-1) ]
    
    if n % 2 == 0: # If an even number of points
        k = np.concatenate( (np.arange(0, n/2+1), np.arange(1-n/2, 0)) )
    else:
        k = np.concatenate( (np.arange(0, (n-1)/2+1), np.arange( -(n-1)/2, 0)) )
    k = 2.*np.pi*k/n
    
    # Modify frequencies by dividing by ik
    f[1:] /= (1j * k[1:]) 
    f[0] = 0. # Set the arbitrary zero-frequency term to zero
    
    return ifft(f).real # Reverse Fourier Transform
   

def pic(f, ncells, L):
    """ f contains the position and velocity of all particles
    """
    nparticles = len(f)/2     # Two values for each particle
    pos = f[0:nparticles] # Position of each particle
    vel = f[nparticles:]      # Velocity of each particle

    dx = L / float(ncells)    # Cell spacing

    # Ensure that pos is between 0 and L
    pos = ((pos % L) + L) % L
    
    # Calculate number density, normalised so 1 when uniform
    density = calc_density(pos, ncells, L)
    
    # Subtract ion density to get total charge density
    rho = density - 1.
    
    # Calculate electric field
    E = -fft_integrate(rho)*dx
    
    # Interpolate E field at particle locations
    accel = -periodic_interp(E, pos/dx)

    # Put back into a single array
    return np.concatenate( (vel, accel) )

####################################################################

def run(pos, vel, L, ncells=None, out=[], output_times=np.linspace(0,20,100), cfl=0.5):
    
    if ncells == None:
        ncells = int(np.sqrt(len(pos))) # A sensible default

    dx = L / float(ncells)
    
    f = np.concatenate( (pos, vel) )   # Starting state
    nparticles = len(pos)
    begintime = t.clock()    
    
    time = 0.0
    for tnext in output_times:
        # Advance to tnext
        stepping = True
        while stepping:
            # Maximum distance a particle can move is one cell
            dt = cfl * dx / max(abs(vel))
            if time + dt >= tnext:
                # Next time will hit or exceed required output time
                stepping = False
                dt = tnext - time
            #print "Time: ", time, dt
            f = rk4step(pic, f, dt, args=(ncells, L))
            time += dt
            
        # Extract position and velocities
        pos = ((f[0:nparticles] % L) + L) % L
        vel = f[nparticles:]
        
        # Send to output functions
        for func in out:
            func(pos, vel, ncells, L, time)
        
    endtime = t.clock()
    print endtime-begintime
    
    return pos, vel

####################################################################
# 
# Output functions and classes
#

class Plot:
    """
    Displays three plots: phase space, charge density, and velocity distribution
    """
    def __init__(self, pos, vel, ncells, L):
        
        d = calc_density(pos, ncells, L)
        vhist, bins  = np.histogram(vel, int(np.sqrt(len(vel))))
        vbins = 0.5*(bins[1:]+bins[:-1])
        
        # Plot initial positions
        if got_gridspec:
            self.fig = plt.figure()
            self.gs = gridspec.GridSpec(4, 4)
            ax = self.fig.add_subplot(self.gs[0:3,0:3])
            self.phase_plot = ax.plot(pos, vel, '.')[0]
            ax.set_title("Phase space")
            
            ax = self.fig.add_subplot(self.gs[3,0:3])
            self.density_plot = ax.plot(np.linspace(0, L, ncells), d)[0]
            
            ax = self.fig.add_subplot(self.gs[0:3,3])
            self.vel_plot = ax.plot(vhist, vbins)[0]
        else:
            self.fig = plt.figure()
            self.phase_plot = plt.plot(pos, vel, '.')[0]
            
            self.fig = plt.figure()
            self.density_plot = plt.plot(np.linspace(0, L, ncells), d)[0]
            
            self.fig = plt.figure()
            self.vel_plot = plt.plot(vhist, vbins)[0]
#        plt.ion()
#        plt.show()
        
    def __call__(self, pos, vel, ncells, L, t):
        d = calc_density(pos, ncells, L)
        vhist, bins  = np.histogram(vel, int(np.sqrt(len(vel))))
        vbins = 0.5*(bins[1:]+bins[:-1])
        
        self.phase_plot.set_data(pos, vel) # Update the plot
        self.density_plot.set_data(np.linspace(0, L, ncells), d)
        self.vel_plot.set_data(vhist, vbins)
#        plt.draw()
#        plt.pause(0.05)
        

class Summary:
    def __init__(self):
        self.t = []
        self.firstharmonic = []
        
    def __call__(self, pos, vel, ncells, L, t):
        # Calculate the charge density
        d = calc_density(pos, ncells, L)
        
        # Amplitude of the first harmonic
        fh = 2.*abs(fft(d)[1]) / float(ncells)
        
        print "Time:", t, "First:", fh
        
        self.t.append(t)
        self.firstharmonic.append(fh)

####################################################################
# 
# Functions to create the initial conditions
#

def landau(npart, L, alpha=0.2):
    """
    Creates the initial conditions for Landau damping
    
    """
    # Start with a uniform distribution of positions
    pos = np.random.uniform(0., L, npart)
    pos0 = pos.copy()
    k = 2.*np.pi / L
    for i in range(10): # Adjust distribution using Newton iterations
        pos -= ( pos + alpha*np.sin(k*pos)/k - pos0 ) / ( 1. + alpha*np.cos(k*pos) )
        
    # Normal velocity distribution
    vel = np.random.normal(0.0, 1.0, npart)
    
    return pos, vel


def twostream(npart, L, vbeam=2):
    # Start with a uniform distribution of positions
    pos = np.random.uniform(0., L, npart)
    # Normal velocity distribution
    vel = np.random.normal(0.0, 1.0, npart)
    
    np2 = int(npart / 2)
    vel[:np2] += vbeam  # Half the particles moving one way
    vel[np2:] -= vbeam  # and half the other
    
    return pos,vel


####################################################################

def setup(mode=1, npart=10000, ncells=20, L = 4.*np.pi):
    # Generate initial condition
    #
    if (mode == 0):
        # 2-stream instability mode
        L = 100
        ncells = 20
        npart = 10000
        pos, vel = twostream(npart, L, 3.)
    elif (mode == 1):
        # Landau damping mode
        L = 4.*np.pi
        ncells = 20
        npart = 20000
        pos, vel = landau(npart, L)
    elif (mode == 2):
        # Landau damping custom mode
        pos, vel = landau(npart, L)
    else:
        pos, vel = landau(npart, L)

    
    # Create some output classes
    p = Plot(pos, vel, ncells, L) # This displays an animated figure
    s = Summary()                 # Calculates, stores and prints summary info
    
    # Run the simulation
    begintime = t.clock()
    pos, vel = run(pos, vel, L, ncells, 
                   out=[p, s],                      # These are called each output
                   output_times=np.linspace(0.,20,100)) # The times to output
    endtime = t.clock()
    
    # Peak finding
    tArr = np.array(s.t)
    fh = np.array(s.firstharmonic)
    peaks = spy.argrelmax(fh)
    
    # Limit to landau daming regime
    lPeaks = peaks[0][0:firstMin(fh[peaks])+1]
    lPeaks = np.insert(lPeaks,0,0)
    
    
    # Summary stores an array of the first-harmonic amplitude
    # Make a semilog plot to see exponential damping
    plt.figure()
    plt.plot(s.t, s.firstharmonic)
    plt.plot(tArr[lPeaks],fh[lPeaks],'x')
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    
    plt.ioff() # This so that the windows stay open
    plt.show()
    return endtime-begintime
    

def repeatedNpartRuns(n=20,Ncell=20,scalefactor=5000):
    data = [[],[]]
    for i in range(n):
        data[0].append(scalefactor*(i+1))
        data[1].append(setup(2,scalefactor*(i+1),Ncell))
    
    plt.figure()
    plt.plot(data[0],data[1])
    plt.xlabel("$N_{part}$")
    plt.ylabel("Simulation Time (seconds)")  
    plt.title("Simulation Time as a Function of $N_{part}$")
    plt.ioff() # This so that the windows stay open
    plt.show()
    
def repeatedNcellRuns(n=20,Npart=10000,scalefactor=5):
    data = [[],[]]
    for i in range(n):
        data[0].append(scalefactor*(i+1))
        data[1].append(setup(2,Npart,scalefactor*(i+1)))
    
    plt.figure()
    plt.plot(data[0],data[1])
    plt.xlabel("$N_{cell}$")
    plt.ylabel("Simulation Time (seconds)")   
    plt.title("Simulation Time as a Function of $N_{cell}$")
    plt.ioff() # This so that the windows stay open
    plt.show()
 
 
####################################################################
   
if __name__ == "__main__":
    setup()
