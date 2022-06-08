"""
hw1_p1.py

Author: Owen Morehead
Date: April 14, 2022

Purpose: AM213B, HW1, Problem 1

"""

# External Modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

# Internal Modules
from Numerical_Schemes import *


#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)


def analytical_sol(t):
    """
    Analytical solution to the two dimensional ODE: dydt = Ay

    A = [[-3, 1],
         [-1, -3]]

    y_0 = [-3,1]^T

    Analytical Solution:
    y(t) = e^(-3t) [sin(t) - 3cos(t)), (cos(t) + 3sin(t))]^T

    Parameters:
    ----------
    t: 1-D array of time values, [N]
    
    Returns:
    -------
    y: 2-D Array of discrete analytical solution values, [N,2]

    """

    N = len(t)           #number of time steps -> grid size

    y = np.zeros((N,2))  #initialize y

    y[:,0] = np.exp(-3*t) * (np.sin(t) - 3*np.cos(t)) #first vector soluiton, y1
    y[:,1] = np.exp(-3*t) * (np.cos(t) + 3*np.sin(t)) #second vector soluiton, y2

    return y


def dydt_p1b(y,t):
    """
    ODE to be solved numerically: dy/dt = Ay

    Parameters
    ----------
    y: 1-D Array of length 2 -- Vector of y values.
    t: Float -- Time coordinate if dealing with time dependent ODE (does not apply here).

    Returns
    -------
    dydt: 1-D Array of length 2
          The derivative solution to the ODE.
          
    """

    dydt = np.dot(A,y)

    return dydt


def AM2_p1b(A, N, dt, t, dydt, y0):
    """
    This function uses the Implicit Adams Moulton method to solve the ODE
    specific to this problem. 

    Parameters
    ----------
    A: 2-D Array, [N,M] -- Matrix of values provided in the problem.
    N: Integer          -- The number of grid points.
    dt: Float           -- Time step.
    t: 1-D Array, [N]   -- Uniformly spaced time grid.
    dydt: Function      -- RHS of the differential equation to be solved.
    y0: 1-D Array, [M]  -- The initial condition to the ODE.
    
    Returns
    -------
    y: 2-D Array, [N,M] -- The solution to the ODE.

    """

    M = len(y0) #number of y vectors in the system

    y = np.zeros((N,M))

    y[0,:] = y0

    for k in range(N-1):

        if k == 0:
        #Using an order 3 explicit method to get y[1,:] (i.e. RK3)
            k1 = dydt(y[k,:],t[k])
            k2 = dydt(y[k,:] + dt*0.5*k1, t[k] + 0.5*dt)
            k3 = dydt(y[k,:] + dt*(-k1 + 2*k2), t[k] + dt)
        
            y[k+1,:] = y[k,:] + dt*((1/6)*k1 + (2/3)*k2 + (1/6)*k3)
        
        else:
        #Use AM2 method for rest of time steps
        #y_k+1 = (1-dt/12 *5*a)^-1 (u_k + dt/12(8*A*u_k - A*u_k-1))

            y[k+1, :] = np.linalg.solve(np.identity(y0.shape[0])-(dt/12)*5*A, y[k,:] + (dt/12)*(8*np.dot(A, y[k,:]) - np.dot(A, y[k-1, :])))
    
    return y


if __name__ == '__main__':


    A = np.array([[-3,1],
                  [-1,-3]])

    T = 10  #total time

    dt = 0.0005  #time step size

    t = np.arange(0,T,dt)  #grid of time values

    N = len(t) #total number of grid points

    y0 = np.array([-3, 1])  #initial condition

    y = analytical_sol(t)

    y_rk3 = RK3(N,dt,t,dydt_p1b,y0)       #RK3 Method

    y_am2 = AM2_p1b(A,N,dt,t,dydt_p1b,y0) #AM2 Method


    #Plotting
    
    #---- Plot of Numerical vs Analytical Solutions ----
    
    fig, ax = plt.subplots()

    ax.plot(t,y,'-o', markersize = 1, label = 'Analytical Soluiton')
    ax.plot(t,y_rk3, '-o', markersize = 1, label = 'RK3 Soluiton')
    ax.plot(t,y_am2, '-o',markersize = 1, label = 'AM2 Solution')

    ax.set_ylabel('y(t)',fontsize=16)
    ax.set_xlabel('t',fontsize=16)
    ax.set_title(f'Solution Comparison: dt = {dt}',fontsize=18)
    ax.legend()
    ax.grid()

    fig.show()
    
    #---- Plot of Error between Numerical and Analytical Solutions ----
    #---- for various values of dt
    
    #dt_vals = [0.5,0.1,0.05,0.005,0.0005,0.00005]
    dt_vals = [0.5,0.1,0.05]
    
    
    fig, axes = plt.subplots(nrows=1,ncols=3)
    axs = axes.flatten()
    
    for j, dt in enumerate(dt_vals):  #Loop over values of dt to calculate errors
    
        t = np.arange(0,T,dt)
        N = len(t)
        
        y0 = np.array([-3,1])
        
        y = analytical_sol(t)
        
        y_rk3 = RK3(N,dt,t,dydt_p1b,y0)
        y_am2 = AM2_p1b(A,N,dt,t,dydt_p1b,y0)
        
        error_rk3 = np.zeros(N)
        error_am2 = np.zeros(N)
        
        for i in range(N):
            error_rk3[i] = np.linalg.norm(y[i,:] - y_rk3[i,:], ord=2)
            error_am2[i] = np.linalg.norm(y[i,:] - y_am2[i,:], ord=2)
            #error_rk3[i] = np.sum((y[i,:]-y_rk3[i,:])**2)
            #error_am2[i] = np.sum((y[i,:]-y_am2[i,:])**2)
            
            
        axs[j].plot(t, error_rk3, 'r', label=f'RK3 Error')
        axs[j].plot(t, error_am2, 'b', label=f'AM2 Error')
  
        axs[j].set_ylabel(f'$||u(t_k) - y(t_k)||_2$',fontsize=15)
        axs[j].set_xlabel('t',fontsize=15)
        axs[j].set_title(f'dt = {dt}',fontsize=15)
        axs[j].legend(title='Method')
        axs[j].grid()
        
        #axs[j].set_yscale('log')
        #axs[j].set_xlim(.5,T)
        #axs[j].set_ylim(1e-50,1)
    
        fig.subplots_adjust(wspace=0.5, hspace=0.5) 
        plt.suptitle("2-Norm Error vs. Time",fontsize=16)
        fig.show()      
        
    
