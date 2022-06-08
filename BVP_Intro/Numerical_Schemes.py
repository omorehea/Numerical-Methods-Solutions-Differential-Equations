"""
Numerical_Schemes.py

Author: Owen Morehead
Date: April 14, 2022

Purpose: Contains a variety of Numerical Schemes to solve a given ODE system.

"""

#External Modules
import numpy as np


def RK3(N, dt, t, dydt, y0):
    """
    This function implements an explicit 3 stage Runge Kutta method to solve a given ODE.
    This algorithm assumes a uniformly spaced grid.
    
    Parameters 
    ----------
    N: Integer          -- The number of grid points.
    dt: Float           -- Time step.
    t: 1-D Array, [N]   -- Uniformly spaced time grid.
    dydt: Function      -- RHS of the differential equation to be solved.
    y0: 1-D Array, [M]  -- The initial condition to the ODE.

    Returns
    -------
    y: 2-D Array, [N,M] -- The discrete solution to the ODE.

    """

    # The length of y0 describes the total number of vectors in the problem
    M = len(y0)
    
    y = np.zeros((N,M))

    y[0,:] = y0  #initial condition

    #loop over time values
    for k in range(N-1):

        #three steps to RK3 method
        k1 = dydt(y[k,:],t[k])
        k2 = dydt(y[k,:] + dt*0.5*k1, t[k] + 0.5*dt)
        k3 = dydt(y[k,:] + dt*(-k1 + 2*k2), t[k] + dt)

        #calcualte next step in y using k1,k2,k3
        y[k+1,:] = y[k,:] + dt*((1/6)*k1 + (2/3)*k2 + (1/6)*k3)

    return y


def Euler_Forward(N,dt,t,dydt,y0):
    """
    This funciton implements the Euler Forward scheme to solve the given ODE.
    This algorithm assumes a uniformly spaced grid.
  
    Parameters 
    ----------
    N: Integer          -- The number of grid points.
    dt: Float           -- Time step.
    t: 1-D Array, [N]   -- Uniformly spaced time grid.
    dydt: Function      -- RHS of the differential equation to be solved.
    y0: 1-D Array, [M]  -- The initial condition to the ODE.

    Returns
    -------
    y: 2-D Array, [N,M] -- The discrete solution to the ODE.

    """ 

    M = len(y0)

    y = np.zeros((N,M))

    y[0,:] = y0

    #calculate next step in y using previous step
    for k in range(N-1):
        y[k+1,:] = y[k,:] + dt*dydt(y[k,:], t[k])

    return y


def Heun_RK2(N,dt,t,dydt,y0):
    """
    This funciton implements the Heun method (explicit RK2) to solve the given ODE.
    This algorithm assumes a uniformly spaced grid.
    
    Parameters 
    ----------
    N: Integer          -- The number of grid points.
    dt: Float           -- Time step.
    t: 1-D Array, [N]   -- Uniformly spaced time grid.
    dydt: Function      -- RHS of the differential equation to be solved.
    y0: 1-D Array, [M]  -- The initial condition to the ODE.

    Returns
    -------
    y: 2-D Array, [N,M] -- The discrete solution to the ODE.

    """    
    
    # The length of y0 describes the total number of vectors to the problem
    M = len(y0)
    
    y = np.zeros((N,M))

    y[0,:] = y0  #initial condition  
  
    #loop over time values
    for k in range(N-1):

        #two steps of RK2
        k1 = dydt(y[k,:],t[k])
        k2 = dydt(y[k,:] + dt*k1, t[k]+dt)
        
        #calcualte next step in y using k1,k2
        y[k+1,:] = y[k,:] + 1/2*dt*(k1 + k2)


    return y


