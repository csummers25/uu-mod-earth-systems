#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualisation routines for plotting code output

"""

from output.visualisation import plotSummary, plotMarkers_lithology, plotTemperature


def makePlots(grid, markers, params, ntstp, t_curr):
    """
    Wrapper function which calls all plotting routines, to simplify calling in 
    the run script.

    Parameters
    ----------
    grid : grid Object
        grid object containing the all the simulation variables on the grid.
    markers : Markers object
        markers object containing the current marker quantities.
    params : Parameters object
        Object containing parameters for the simulation.
    ntstp: INT
        Current timestep number
    t_curr : FLOAT
        Current simulation time.

    Returns
    -------
    None.

    """
    xlims = (0,params.xsize)
    ylims = (params.ysize,0)
    title = 'Time: %.3f Myr'%(t_curr*1e-6/(365.25*24*3600))
    
    # resolution for the markers plots
    xres = 401
    
    plotSummary(grid, params, ntstp, t_curr, xlims, ylims, title, rhomin=3200, rhomax=3300, vmin=19, vmax=21)
    plotMarkers_lithology(params, markers, grid, ntstp, t_curr, xlims, ylims, title, xres, height=6)
    plotTemperature(grid, params, ntstp, t_curr, xlims, ylims, title, height=6, Tmin=1000, Tmax=1100)


    
    
    
