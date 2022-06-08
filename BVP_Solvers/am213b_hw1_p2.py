"""
hw1_p2.py

Author: Owen Morehead
Date: April 13, 2022

Purpose: AM213B, HW1, Problem 2

Computes the numerical solution to the ODE energy-preserving system.
Plots the solutions in x and theta, and the Hamiltonian for Numerical schemes:
-- Euler Forward
-- Heun Method (explicit RK2)
-- Implicit Midpoint Method

"""

# External Modules
import numpy as np
import matplotlib.pyplot as plt

# Internal Modules
from Numerical_Schemes import *


def dydt_p2(y,t):
    """
    System of differential equations of the Hamiltonial system to be solved.
    Hamiltonian System in vector components:
    y1, y2, y3, y4 = [x, theta, p_x, p_theta]
        
        dy_0/dt = y_2/m
        dy_1/dt = y_3/(m(l_0 + y_0)**2)
        dy_2/dt = y_3/(m(l_0 + y_0)**3) + (mg)cos(y_1) - (kappa)y_0
        dy_3/dt = -(mg)(l_0 + y_0)sin(y_1)


    Parameters
    ----------
    y: 1-D Array of length 4    -- Vector of y values.
    t: Float -- Time coordinate if dealing with time dependent ODE (does not apply here).

    Returns
    -------
    dydt: 1-D Array of length 4 -- Hamiltonians equations of motion.

    """

    dydt = np.zeros(len(y))

    dydt[0] = y[2]/m
    dydt[1] = y[3]/(m*(l_0 + y[0])**2)
    dydt[2] = y[3]**2/(m*(l_0 + y[0])**3) + m*g*np.cos(y[1]) - kappa*y[0]
    dydt[3] = -m*g*(l_0 + y[0])*np.sin(y[1])
    
    return dydt

def Hamiltonian(y):
    """
    Hamiltonian function of the provided dynamical system.
    
    Parameters
    ----------
    y: 2-D Array with 4 columns, [len(t),4] -- Array of y values.
     
    Returns
    -------
    H: 1-D Array, [len(t)] -- Hamiltonian evaluated at timestep t_k.
    
    """
    
    H = y[:,2]**2/(2*m) + y[:,3]**2/(2*m*(l_0+y[:,0])**2) -m*g*(l_0+y[:,0])*np.cos(y[:,1]) + 1/2*kappa*y[:,0]**2
    
    return H

def Jacobian(y,dt):
    """
    Jacobian matrix for the ODE system.
    
    Parameters
    ----------
    y: 1-D Array of length 4, [4] -- Array of y values at time t_k.
    dt: Float                     -- value of timestep.
     
    Returns
    -------
    jac: 4 x 4 Array, [len(t)]    -- Jacobian matrix evaluated at timestep t_k.
    
    """    
    

    jac = [[0, 0, 1/m, 0], 
           [-(2*y[3]) / (m*(l_0+y[0])**3), 0, 0, 1/(m*(l_0+y[0])**2)], 
           [-(3*y[3]**2)/(m*(l_0+y[0])**4) - kappa, -m*g*np.sin(y[1]), 0, (2*y[3])/(m*(l_0 + y[0])**3)], 
           [-m*g*np.sin(y[1]), -m*g*(l_0+y[0])*np.cos(y[1]), 0, 0]]

    jac = 1/2*dt*np.array(jac) #multiply by 0.5dt from scheme equation

    return jac


def Midpoint_Newton(N,dt,t,dydt,y0,tol=10e-5,itermax=50):
    """
    This algorithm uses the Newton's iterative method to implement the 
    Implicit Midpoint Method to solve the Hamiltonian system provided.
    This problem uses the analytically derived Jacobian.
    
    This algorithm assumes a uniformly spaced grid.

    Parameters
    ----------
    N: Integer          -- The number of grid points.
    dt: Float           -- Time step.
    t: 1-D Array, [N]   -- Uniformly spaced time grid.
    dydt: Function      -- RHS of the differential equation to be solved.
    y0: 1-D Array, [M]  -- The initial condition to the ODE.
    tol: Float, optional parameter, default = 10e-5.
    itermax: Float, optional parameter, default = 50 iterations.
    
    Returns
    -------
    y: 2-D Array, [N,M] -- The discrete solution to the ODE.
    
    """
    
    M = len(y0)          #number of solution vectors
    y = np.zeros((N,M))  
    
    I = np.identity(4)   #4 x 4 identity matrix
    
    y[0,:] = y0          #Initial condition
    
    #loop over time values to calculate discrete solution
    for k in range(N-1):
        #implement Newton's method on current iteration step
        
        test = 1  #test placeholder value for convergence
        iters = 0 #hold the iteration count
        
        #while loop continues until u_k+1 is below tolerance
        while (test>tol):
            y_prev = y[k+1,:] #store previous y
            
            if (iters>itermax):
                print("Max iteration count reached")
                exit()
                
            if iters == 0:        #initial guess to k+1 step
                y[k+1,:] = y[k,:]
    
            else: #solve the nonlinear system using Newton's method to calculate y[k+1,:]
                x = np.linalg.solve(I-Jacobian(y[k+1,:],dt),y[k,:] + dt*dydt(1/2*(y[k,:] + y[k+1,:]),t[k]+1/2*dt) - y[k+1,:])
                    
                y[k+1,:] = x + y[k+1,:]
                
                #calculate the norm between old and new iteration values 
                # to check convergence
                test = np.linalg.norm(y[k+1,:]-y_prev, ord=2)
    
            iters += 1
            
    return y
    

if __name__ == '__main__':

    #--- Global Parameters ---

    g = 9.8     # m/s^2
    l_0 = 1     # m
    kappa = 10  # kg
    m = 0.5     # N/m


    y0 = np.zeros(4)
    
    #--- Initial Conditions ---
    y0[0] = 0            #x(0) = 0
    y0[1] = (3/4)*np.pi  #theta(0) = 3pi/4
    y0[2] = 0            #p_x(0) = 0
    y0[3] = 0            #p_theta(0) = 0


    T = 60

    dt = 4*10**(-2)

    t = np.arange(0,T,dt)   #array of uniformly spaced time values
   
    N = len(t)              #size of time array -> total number of discrete solution points

    #--- Solving using Euler Forward, Heun Method, and Implicit Midpoint ---
    y_EF = Euler_Forward(N,dt,t,dydt_p2,y0)
    y_H = Heun_RK2(N,dt,t,dydt_p2,y0)
    y_IM = Midpoint_Newton(N,dt,t,dydt_p2,y0)
    
    #--- Hamiltonian of each solution scheme ---
    Ham_EF = Hamiltonian(y_EF)
    Ham_H = Hamiltonian(y_H)
    Ham_IM = Hamiltonian(y_IM)


    #--- Plotting ---
    
    # Euler-Forward Plots ------------------

    fig, axs = plt.subplots(nrows=1,ncols=3)
    
    ax = axs[0] #plotting x and theta vs t for Euler Method
    
    ax.plot(t,y_EF[:,0],"r",label='x')
    ax.plot(t,y_EF[:,1],"b",label=r'$\theta$')
    
    ax.set_title("Euler Forward Method")
    ax.set_xlabel("t (s)",fontsize=14)
    ax.legend()
    ax.grid()

    
    # Heun Method Plots ------------------

    ax = axs[1] #plotting x and theta vs t for Heun Method
    ax.plot(t,y_H[:,0],"r",label='x')
    ax.plot(t,y_H[:,1],"b",label=r'$\theta$')
    
    ax.set_title("Heun Method")
    ax.set_xlabel("t (s)",fontsize=14)
    ax.legend()
    ax.grid()

    # Implicit Midpoint Method Plots ------------------

    ax = axs[2] #plotting x and theta vs t for Implicit Midpoint Method
    ax.plot(t,y_IM[:,0],"r",label='x')
    ax.plot(t,y_IM[:,1],"b",label=r'$\theta$')

    ax.set_title("Implicit Midpoint Method")
    ax.set_xlabel("t (s)",fontsize=14)
    ax.legend()
    ax.grid()

    plt.suptitle(fr"x and $\theta$ vs. Time -- $\Delta$t = {dt}",fontsize=16)
    fig.subplots_adjust(wspace=0.4, hspace=0.4) 
    fig.show()
    
    
    # Hamiltonian Plot for all Methods
    # ------------------
    
    fig, ax = plt.subplots()
    ax.plot(t,Ham_EF,label='Euler Forward')
    ax.plot(t,Ham_H,label='Heun Method')
    ax.plot(t,Ham_IM, label='Implicit Midpoint')
    
    ax.set_title(f"Hamiltonian Funciton vs. Time -- $\Delta$t = {dt}",fontsize=16)
    ax.set_xlabel("t (s)",fontsize=14)
    ax.legend()
    ax.grid()
    #ax.set_ylim(0,500)


    fig.show()
    
    
    
