# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:12:34 2014

kgsimulate simulates a solution for the 2nd order Klein-Gordon wave equation 
using a total spectral method. 

The simulation is ran in k-space (fourier transform of real space), in this 
form the (spectral) solution completely diaganolizes and we get Nk^3 
uncoupled equations for the fourier coefficients of our Electric Field.

Because of this, time-stepping is quite trivial, and the calculation is very
quick (but speed falls off in third order, which is usual of 3d problems)

Current problems to be addressed:
    Boundaries! With a discrete # of k values our problem is periodic, i.e. the
    simulation is running in a "box" and the waves will reflect on the edges.

@author: Stephen Carr, Brian Dawes
"""

import numpy as np
from numpy import pi

"""
Here we define our global variables
"""

#timestep h
h = .1

#total amount of timesteps
tsteps = 250

#number of mesh points in k-space
#note, odd number means we have k = 0
Nk = 51

#mass
m = 1

"""
Next two methods are for our source, f
"""

#fourier coefficients of our source in the d'th component (d = 3 is z, etc)
def f_k(t,kx,ky,kz,d):
    #if t == 0:
        if d == 1:
            if t <= 20*h:
                return -2*1j*kx/(1+np.exp(-0.25*(t/h-10)))*np.exp(-kx*0.2*pi*1j*cos(2*pi*0.1*t)-ky*0.2*pi*1j*sin(2*pi*0.1*t))
            return -2*1j*kx*np.exp(-kx*0.2*pi*1j*cos(2*pi*0.1*t)-ky*0.2*pi*1j*sin(2*pi*0.1*t))
            
        if d == 2:
            if t <= 20*h:
                return -2*1j*ky/(1+np.exp(-0.25*(t/h-10)))*np.exp(-kx*0.2*pi*1j*cos(2*pi*0.1*t)-ky*0.2*pi*1j*sin(2*pi*0.1*t))
            return -2*1j*ky*np.exp(-kx*0.2*pi*1j*cos(2*pi*0.1*t)-ky*0.2*pi*1j*sin(2*pi*0.1*t))
            
        if d == 3:  
            if t <= 20*h:
               return -2*1j*kz/(1+np.exp(-0.25*(t/h-10)))
            return -2*1j*kz
    #else:
    #    return 0

#matrix of the fourier coefficients in k space
def f(t,d):
    
    f_arr = np.zeros((Nk,Nk,Nk),dtype = complex)
    
    # k's are ordered from 0 to 2pi then negatives to -2pi
    # this is to match fft convention
    for i in xrange(0,Nk):
        kx = 4*pi*i/(Nk-1)
        if kx>2*pi:
            kx = kx-4*pi
        for j in xrange(0,Nk):
            ky = 4*pi*j/(Nk-1)
            if ky>2*pi:
                ky = ky-4*pi
            for l in xrange(0,Nk):
                kz = 4*pi*l/(Nk-1)
                if kz>2*pi:
                    kz = kz-4*pi
                f_arr[kx,ky,kz] = f_k(t,kx,ky,kz,d)
                
    return f_arr

'''
knock knock
whos there?
RAVIOLI
RAVIOLI
RAVIOLI
'''

"""
spectral(d) runs our spectral simulation for the d'th componenet of the E field
"""

def spectral(d):
    
    """ 
    define variables for our time stepping process
    for now assume inital conditions are zero
    Also set up k and mass arrays for the time step
    """
    
    #initial conditions (phi_last should always = phi_curr I think)
    phi_curr = np.zeros((Nk,Nk,Nk),dtype = complex)
    phi_last = phi_curr
    
    #k^2 term as an array
    k_arr = np.zeros((Nk,Nk,Nk),dtype = complex)
    # k's are ordered from 0 to 2pi then negatives to -2pi
    # this is to match fft convention
    for i in xrange(0,Nk):
        kx = 4*pi*i/(Nk-1)
        if kx>2*pi:
            kx = kx-4*pi
        for j in xrange(0,Nk):
            ky = 4*pi*j/(Nk-1)
            if ky>2*pi:
                ky = ky-4*pi
            for l in xrange(0,Nk):
                kz = 4*pi*l/(Nk-1)
                if kz>2*pi:
                    kz = kz-4*pi
                k_arr[kx,ky,kz] = kx**2 + ky**2 + kz**2
    
    solution = np.zeros((tsteps+1,Nk,Nk,Nk),dtype = complex)
    
    #running the simulation, only 5 lines!
    for t in range(0,tsteps):            
        phi_next = 2*phi_curr - phi_last + (h**2)*(f(t*h,d) - (k_arr+m**2)*phi_curr)
        solution[t,:,:,:] = phi_next
        phi_last = phi_curr
        phi_curr = phi_next
        print "Time Step " + str(t) + " done."
    
    return solution
    
def main():
    
    #Run the simulation in the d = 1,2,3 cases (corresponds to component x,y,z)
    print 'Running K_x'    
    K_x = spectral(1)
    print 'Running K_y'  
    K_y = spectral(2)
    print 'Running K_z'  
    K_z = spectral(3)
    
    #save our results
    np.save("K_x", K_x)
    np.save("K_y",K_y)
    np.save("K_z",K_z)
    
    #save the parameters (makes plotting easier)
    params = np.array([h, tsteps, Nk, m])
    np.save("params",params)
    
    print "Done!"
    
if __name__ == "__main__":
    main()