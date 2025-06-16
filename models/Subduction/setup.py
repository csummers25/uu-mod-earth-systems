#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File which contains the setup for a lithospheric extension model

"""

import numpy as np
import pathlib
import math

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
    
    for j in range(0,markers.xnum):
        for i in range(0,markers.ynum):
            
            # set coordinates as grid + small random displacement
            markers.x[mm] = (j + np.random.random())*mxstp
            markers.y[mm] = (i + np.random.random())*mystp
            
            # now define rock type based on location, to give our initial setup
            
            # asthenosphere
            markers.id[mm] = 5
            
            # sticky air
            if (markers.y[mm] <= 10000):
                markers.id[mm] = 0

            # Sediments
            if (markers.y[mm] > 10000 and markers.y[mm] <= 11000):
                markers.id[mm] = 1
            # Basaltic crust 
            if (markers.y[mm] > 11000 and markers.y[mm] <= 13000):
                markers.id[mm] = 2
            # Gabbroic crust
            if (markers.y[mm] > 13000 and markers.y[mm] <= 18000):
                markers.id[mm] = 3
            # Lithospheric mantle
            if (markers.y[mm] > 18000 and markers.y[mm] < 90000):
                markers.id[mm] = 4
            # Continent
            # Upper continental crust
            if (markers.x[mm] < 200000 and markers.y[mm] > 7000 and markers.y[mm] < 27000):
                markers.id[mm] = 7
            if (markers.x[mm] >= 200000 and markers.x[mm] < 250000 and  markers.y[mm] > 7000 + (markers.x[mm]-200000)/50000*3000 and markers.y[mm] < 32000 - (markers.x[mm]-200000)/50000*14000):
                markers.id[mm] = 7
            # Used for visualization of the model
            if (markers.x[mm] < 200000 and markers.y[mm] > 11000 and markers.y[mm] < 15000):
                markers.id[mm] = 8
            if (markers.x[mm] < 200000 and markers.y[mm] > 19000 and markers.y[mm] < 23000):
                markers.id[mm] = 8
            # lower cont crust
            if (markers.x[mm] < 200000 and markers.y[mm] > 27000 and markers.y[mm] < 42000):
                markers.id[mm] = 9
            if (markers.x[mm] < 200000 and markers.y[mm] > 32000 and markers.y[mm] < 37000):
                markers.id[mm] = 10
            if (markers.x[mm] >= 200000 and markers.x[mm] < 250000 and markers.y[mm] > 27000 - (markers.x[mm]-200000)/50000*14000 and markers.y[mm] < 42000 - (markers.x[mm]-200000)/50000*24000):
                markers.id[mm] = 9
    	    # Weak zone in the mantle
            if (markers.y[mm] > 28000 and markers.y[mm] < 90000 and markers.x[mm] > 240000 - (markers.y[mm]-18000)/24000*50000 and markers.x[mm] < 260000 - (markers.y[mm]-18000)/24000*50000):
                markers.id[mm] = 6            
	         # Weak zone in the oceanic crust
            if (markers.y[mm] > 20000 and markers.y[mm] <= 18000 and markers.x[mm] > 240000 - (markers.y[mm]-18000)/24000*50000 and markers.x[mm] < 260000 - (markers.y[mm]-18000)/24000*155000):
                markers.id[mm] = 2
	        # Asthenosphere below the oceanic plate
            if (markers.y[mm] > 70000 and markers.y[mm] <= 100000 and markers.x[mm] > 410000 - (markers.y[mm]-8000)/24000*75000):
                markers.id[mm] = 5


            # lithospheric mantle, default is 5, we add stripes of 4
            #if (markers.y[mm]>52000 and markers.y[mm]<= 57000):
            #    markers.id[mm] = 4
            #if (markers.y[mm]>62000 and markers.y[mm]<= 67000):
            #      markers.id[mm] = 4
            #if (markers.y[mm]>72000 and markers.y[mm]<= 77000):
            #    markers.id[mm] = 4
            #if (markers.y[mm]>82000 and markers.y[mm]<= 87000):
            #    markers.id[mm] = 4
            #if (markers.y[mm]>92000 and markers.y[mm]<= 97000):
            #    markers.id[mm] = 4
            
            # initial temperature structure
            # adiabatic T gradient in asthenosphere
            dtdy = 0.5/1000 # K/m
            markers.T[mm] = params.T_bot - dtdy*(ysize - markers.y[mm])
            
            # if in the air
            if (markers.id[mm]==0):
                # cons T
                markers.T[mm] = params.T_top
            
	        # Oceanic geotherm
            age = 3e7*(365.25*24*3600)      #Oceanic plate age, s
            y_asth = 92000   	           #Bottom of the lithosphere
            T_asth = params.T_bot - dtdy*(ysize - y_asth)  # T of astenosphere at y=yast
            kappa = 1e-6		     # Thermal diffusivity of the mantle, m^2/s
    
	        # After Turcotte & Schubert (2002) 
            if (markers.x[mm] > 250000 and markers.y[mm] > 10000 and markers.y[mm] < y_asth):  # markers.y was creater than 10000
                # T difference at the bottom of the oceanic plate
                dt = -(params.T_top - T_asth)*(1-math.erf((y_asth-10000)/2/(kappa*age)**0.5))
                markers.T[mm] = T_asth+dt+(params.T_top-T_asth-dt)*(1-math.erf((markers.y[mm]-10000)/2/(kappa*age)**0.5))

            # linear continental geotherm
            if (markers.x[mm] <= 200000 and markers.y[mm] > 7000 and markers.y[mm]<y_asth):
                markers.T[mm] = params.T_top + (T_asth - params.T_top)*(markers.y[mm] - 7000)/(y_asth - 7000)
            
	        # Transient left continent --> ocean geotherm
            if (markers.x[mm] > 200000 and markers.x[mm] < 250000 and markers.y[mm] > 7000 + (markers.x[mm]-200000)/50000*3000 and markers.y[mm] < y_asth):
		        # Continental geotherm
                T_cont = params.T_top + (T_asth-params.T_top)*(markers.y[mm]-(7000+(markers.x[mm]-200000)/50000*3000))/(y_asth-7000)
		        # T difference at the bottom of the oceanic and continental plates
                dt = -(params.T_top - T_asth)*(1-math.erf((y_asth-10000)/2/(kappa*age)**0.5))
		        # Oceanic geotherm
                T_ocea = T_asth+dt+(params.T_top-T_asth-dt)*(1-math.erf((markers.y[mm]-(7000+(markers.x[mm]-200000)/50000*3000))/2/(kappa*age)**0.5))
		        # Linear lateral transition
                mwt = (markers.x[mm]-200000)/50000
		        # Transitional temperate
                markers.T[mm] = T_cont*(1-mwt) + T_ocea*mwt

            # update marker index
            mm +=1
    
    
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
        B_intern[3] = prescribed x-velocity value.
        B_intern[4] = y-index of vy nodes with prescribed velocity (-1 is not in use) 
        B_intern[5-6] = min/max x-index of the wall. 
        B_intern[7] = prescribed y-velocity value.
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

    # additional model options 
    # initial system size
    xsize0 = 760000      #600000    
    ysize0 = 200000
    
    xsize = xsize0
    ysize = ysize0

    # set resolution
    xnum = 104+4
    ynum = 35


    # instantiate/load material properties object
    # file path must be from top directory (as that is where the fn is called from!)
    matData = np.loadtxt('./models/Subduction/material_properties.txt', skiprows=3, delimiter=",")
    materials = Materials(matData)

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
    B_bottom[:,1] = 1
    B_bottom[:,2] = 0

    B_left = np.zeros((ynum+1,4))
    B_left[:,0] = 0 
    B_left[:,3] = 1

    B_right = np.zeros((ynum+1,4))
    B_right[:,0] = 0 
    B_right[:,3] = 1

    # optional internal boundary, switched off
    B_intern = np.zeros(8)
    B_intern[0] = 99  
    B_intern[1] = 20
    B_intern[2] = 23
    B_intern[3] = (-7.5*1e-2)/(365.25*24*3600)  # convert to m/s
    B_intern[4] = -1 
    B_intern[5] = 0
    B_intern[6] = 0
    B_intern[7] = 0  

    # temperature BCs
    BT_top = np.zeros((xnum, 2))
    BT_bottom = np.zeros((xnum, 2))
    BT_left = np.zeros((ynum, 2))
    BT_right = np.zeros((ynum, 2))

    # upper and lower  = fixed T
    BT_top[:,0] = params.T_top
    BT_bottom[:,0] = params.T_bot

    # left and right = insulating?
    BT_left[:,1] = 1
    BT_right[:,1] = 1

    ###########################################################################
    # create grid object
    grid = Grid(xnum, ynum)

    # define grid points for (potentially) unevenly spaced grid
    gridSpacings(params.bx, params.by, params.Nx, params.Ny, params.non_uni_xsize, xsize, ysize, grid, 0)

    ############################################################################
    # create markers object
    mnumx = 600
    mnumy = 200
    markers = Markers(mnumx, mnumy)

    # initialize markers
    initialize_markers(markers, materials, params, xsize, ysize)
    
    return params, grid, materials, markers, xsize, ysize, P_first, B_top, B_bottom,\
           B_left, B_right, B_intern, BT_top, BT_bottom, BT_left, BT_right
    
    
    











    
