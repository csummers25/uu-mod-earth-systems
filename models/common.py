#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

common functions that can be used by many models

"""

def uniformGrid(grid, xsize, ysize):
    '''
    Calculates the grid node positions and spacings.
    
    This default version implements a fixed, uniform grid, with x=0, y=0 as the starting point.

    Parameters
    ----------
    grid : OBJ
        The grid object into which the new node positions will be written.
    xsize : FLOAT
        The physical size of the simulation domain in the x direction.
    ysize : FLOAT
        The physical size of the simulation domain in the y direction.   
    
    Returns
    -------
    None.

    '''
    
    
    xnum = grid.xnum
    ynum = grid.ynum
    
    dx = xsize/(xnum-1)
    dy = ysize/(ynum-1)
    
    # Simple, uniform grid
    # horizontal grid
    grid.x[0] = 0
    for i in range(1,xnum):
        grid.x[i] = grid.x[i-1] + dx
        
    # vertical grid
    grid.y[0] = 0
    for i in range(1,ynum):
        grid.y[i] = grid.y[i-1] + dy

    # set the array of grid spacings with the new values
    grid.set_spacings()

    # also set the centered node positions + spacings
    grid.set_centered_nodes()