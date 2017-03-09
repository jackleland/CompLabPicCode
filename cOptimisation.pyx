# -*- coding: utf-8 -*-

def calc_density(position, ncells, L):
    """ Optimised method for calculating charge density given 
    particle positions.
    
    INPUT
      position  - Array of positions, one for each particle
                  assumed to be between 0 and L
      ncells    - Number of cells
      L         - Length of the domain

    OUTPUT
      density   - contains 1 if evenly distributed
    """
    
    density = [0]*ncells
        
    nparticles = len(position)
    
    dx = L / ncells       # Uniform cell spacing
    pos = position/ dx
    
    # original algorithm
    for p in position:    # Loop over all the particles, converting position into a cell number
        pos = p / dx        
        plower = int(pos)        # Cell to the left (rounding down)
        offset = pos - plower    # Offset from the left
        density[plower] += 1. - offset
        density[(plower + 1) % ncells] += offset
    
    # nparticles now distributed amongst ncells
    for i in range(ncells):
        density[i] *= float(ncells) / float(nparticles)  # Make average density equal to 1
        
    return density