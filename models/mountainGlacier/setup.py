#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File which contains the setup for a lithospheric extension model

"""

import numpy as np
import pathlib

from dataStructures import Markers, Parameters, Grid, Materials
from physics.grid_fns import gridSpacings


def initialize_markers(markers, materials, params, xsize, ysize):
    '''
    Initialize the positions, material ID and temperature of the markers.

    Parameters
    ----------
    markers : Markers Object
        An empty markers object, which will be filled by this function.
    materials : Materials object
        Materials object filled from file with the required materials.
    params : Parameters Object
        Contains the simulation parameters
    xsize : FLOAT
        Physical x-size of the simulation domain.
    ysize : FLOAT
        Physical y-size of the simulation domain.

    Returns
    -------
    None.

    '''
    
    np.random.seed(1337)
    # set the rough grid spacing for the markers
    mxstp = xsize/markers.xnum
    mystp = ysize/markers.ynum
    
    # marker number index
    mm = 0
    xx, yy = np.meshgrid(markers.xnum, markers.ynum)
    points = np.vstack((xx.ravel(), yy.ravel())).T
    
    for j in range(0,markers.xnum):
        for i in range(0,markers.ynum):
            
            # set coordinates as grid + small random displacement
            markers.x[mm] = (j + np.random.random())*mxstp
            markers.y[mm] = (i + np.random.random())*mystp
            

            

            # bedrock
            if markers.y[mm] >= topography_curve(markers.x[mm], xsize, ysize):
                markers.id[mm] = 2
                markers.T[mm] = 273

            # glacier
            elif markers.y[mm] >= ice_curve(markers.x[mm], xsize, ysize):
                markers.id[mm] = 1
                markers.T[mm] = 263
            
            # air but lapse rate doesnt really seem to work how I want it to work
            else: 
                dtdy = 0.65/1000 # approximate adiabetic lapse rate for the air K/m
                markers.id[mm] = 0 
                markers.T[mm] = 253 - dtdy*(ysize - markers.y[mm])

        


            # grid somehow stops at 1/3 and air temp is not getting in correctly

            # update marker index
            mm +=1



def topography_curve(x, xsize, ysize): 
    '''
    Sets up the curve of the bedrock on which the glacier lays, takes form of a concave up decreasing curve
    Hard-coded, dependent on grid size but have not yet implemented this
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    a = ysize / (xsize ** 2)  # ensures curve drops from ysize to 0
    curve_y =  a * (x)**2  # h = 0, k = 0
    return curve_y 



def ice_curve(x, xsize, ysize):
    '''
    Sets up the curve of the initial height of the glacier, takes the form of a concave down decreasing curve
    Hard-coded, dependent on grid size but have not yet implemented this 
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    a = ysize / (xsize ** 2) 
    curve_y = ysize / (xsize ** 2) * x**2 - ysize/7
    
    return curve_y 

    
def initializeModel():
    '''
    Sets up the initial state of the model, including BCs and output settings.

    Returns
    -------
    params : Parameters object
        Model physical and numerical parameters.
    grid : Grid object
        Initialised Grid object.
    materials : Materials
        Materials object initialsed with required material properties.
    markers : Markers
        Initialized markers object.
    xsize : INT
        x-resolution of model domain.
    ysize : INT
        y-resolution of model domain.
    P_first : ARRAY
        Array with 2 entries, specifying pressure BC.
    B_top : ARRAY
        Boundary conditions at the top of the grid. Array has 4 columns, 
        values in each are defined as: 
        vx[0,j] = B_top[j,0] + vx[1,j]*B_top[j,1]
        vy[0,j] = B_top[j,2] + vy[1,j]*B_top[j,3]
    B_bottom : ARRAY
        Boundary conditions at the bottom of the grid. Array has 4 columns, 
        values in each are defined as: 
        vx[0,j] = B_bot[j,0] + vx[1,j]*B_bot[j,1]
        vy[0,j] = B_bot[j,2] + vy[1,j]*B_bot[j,3]
    B_left : ARRAY
        Boundary conditions at the left of the grid. Array has 4 columns, 
        values in each are defined as: 
        vx[0,i] = B_left[i,0] + vx[1,i]*B_left[i,1]
        vy[0,i] = B_left[j,2] + vy[1,i]*B_left[i,3]
    B_right : ARRAY
        Boundary conditions at the left of the grid. Array has 4 columns, 
        values in each are defined as: 
        vx[0,i] = B_right[i,0] + vx[1,i]*B_right[i,1]
        vy[0,i] = B_right[j,2] + vy[1,i]*B_right[i,3]
    B_intern : ARRAY
        Array defining optional internal boundary eg. moving wall. Format is:
        B_intern[0] = x-index of vx nodes with prescribed velocity (-1 is not in use)
        B_intern[1-2] = min/max y-index of the wall
        B_intern[4] = prescribed x-velocity value.
        B_intern[5] = y-index of vy nodes with prescribed velocity (-1 is not in use)
        B_intern[6-7] = min/max x-index of the wall
        B_intern[8] = prescribed y-velocity value.
    BT_top : ARRAY
        Top temperature BCs.  Array has 2 columns, values in each are defined as:
        T[i,j] = BT_top[0] + BT_top[1]*T[i+1,j]
    BT_bottom : ARRAY
        Bottom temperature BCs.  Array has 2 columns, values in each are defined as:
        T[i,j] = BT_bottom[0] + BT_bottom[1]*T[i-1,j]
    BT_left : ARRAY
        Left temperature BCs.  Array has 2 columns, values in each are defined as:
        T[i,j] = BT_left[0] + BT_left[1]*T[i,j+1]
    BT_right : ARRAY
        Right temperature BCs.  Array has 2 columns, values in each are defined as:
        T[i,j] = BT_right[0] + BT_right[1]*T[i,j-1]

    '''
    
    
    # instantiate a pre-populated parameters object
    params = Parameters()
    
    params.T_top = 273 - 20 # Top of model, temperature is 20 degrees below zero
    params.T_bot = 273 + 2  # Bottom of model, temperature is 2 degrees above zero 
    
    # We must somewhere assume that bedrock has constant temperature and temperature of bedrock remains constant
    
    
    ## time 
    params.t_end = 100 * 365.25 * 3600 * 24  # set end time at t = 100 years for now
    params.tstp_max = 365.25 * 3600 * 24 # max 1 year timestep 
    
    # additional model options 
    # initial system size
    xsize0 = 400000
    ysize0 = 300000
    
    xsize = xsize0
    ysize = ysize0

    # set resolution
    xnum = 610
    ynum = 61


    # instantiate/load material properties object
    # file path must be from top directory (as that is where the fn is called from!)
    matData = np.loadtxt('./models/mountainGlacier/material_properties.txt', skiprows=3, delimiter=",")
    materials = Materials(matData)
    
    
    # For the materials, we in fact only need 3? Air, ice and bedrock 
    

    # create directories for output of figures and data (not atm)
    pathlib.Path('./Figures').mkdir(exist_ok=True)
    pathlib.Path('./Output').mkdir(exist_ok=True)

    ###########################################################################    
    # Boundary conditions
    # pressure BCs
    P_first = np.array([0,1e5])

    # velocity BCs
    B_top = np.zeros((xnum+1,4))
    B_top[:,1] = 1

    B_bottom = np.zeros((xnum+1,4))
    # B_bottom[:,1] = 1
    # B_bottom[:,2] = -params.v_ext/xsize * ysize

    B_left = np.zeros((ynum+1,4))
    # B_left[:,0] = -params.v_ext/2
    # B_left[:,3] = 1

    B_right = np.zeros((ynum+1,4))
    B_right[:,0] = params.v_ext/2
    B_right[:,3] = 1

    # optional internal boundary, switched off
    B_intern = np.zeros(8)
    B_intern[0] = -1
    B_intern[4] = -1
    
    # temperature BCs
    BT_top = np.zeros((xnum, 2))
    BT_bottom = np.zeros((xnum, 2))
    BT_left = np.zeros((ynum, 2))
    BT_right = np.zeros((ynum, 2))

    # upper and lower  = fixed T
    BT_top[:,0] = params.T_top
    BT_bottom[:,0] = params.T_bot

    # left and right = insulating?
    # BT_left[:,1] = 1
    # BT_right[:,1] = 1

    ###########################################################################
    # create grid object
    grid = Grid(xnum, ynum)

    # define grid points for (potentially) unevenly spaced grid
    gridSpacings(params.bx, params.by, params.Nx, params.Ny, params.non_uni_xsize, xsize, ysize, grid, 0)

    ############################################################################
    # create markers object
    mnumx = 400
    mnumy = 300
    markers = Markers(mnumx, mnumy)

    # initialize markers
    initialize_markers(markers, materials, params, xsize, ysize)
    
    return params, grid, materials, markers, xsize, ysize, P_first, B_top, B_bottom,\
           B_left, B_right, B_intern, BT_top, BT_bottom, BT_left, BT_right
    
    
    











    