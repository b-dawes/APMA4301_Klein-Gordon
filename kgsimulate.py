# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:12:34 2014

@author: tfuser
"""
import numpy as np

#timestep h
h = 0.05

#total amount of timesteps
tsteps = 20

#mesh spacing in k-space
dk = 1

#number of mesh points in k-space
#note, odd number means we have k = 0
Nk = 101

k_max = dk*((Nk-1)/2)

#mass
m = 1


def f_k(t,x,y,z,d):
    if d == 1:
        return -4*pi*x
        
    if d == 2:
        return -4*pi*y
        
    if d == 3:  
        return -4*pi*z

def f(t,d):
    f_arr = np.zeros((Nk,Nk,Nk))
    
    for x in range(-k_max,k_max,dk):
        for y in range(-k_max,k_max,dk):
            for z in range(-k_max,k_max,dk):
                f_arr[x,y,z] = f_k(t,x,y,z,d)
                
    return f_arr


def spectral(d):
    
    """ define variables for our time stepping process
        for now assume inital conditions are zero
        Also set up k and mass arrays for the time step
    """
    
    phi_curr = np.zeros((Nk,Nk,Nk))
    phi_last = np.zeros((Nk,Nk,Nk))
    
    k_arr = np.zeros((Nk,Nk,Nk))
    for x in range(-k_max,k_max,dk):
        for y in range(-k_max,k_max,dk):
            for z in range(-k_max,k_max,dk):
                k_arr[x,y,z] = x**2 + y**2 + z**2
    
    solution = np.zeros((tsteps+1,Nk,Nk,Nk))
    
    for t in range(0,tsteps):            
        phi_next = 2*phi_curr - phi_last + (h**2)*(f(t*h,d) - (k_arr+m**2)*phi_curr)
        solution[t,:,:,:] = phi_next
        phi_last = phi_curr
        phi_curr = phi_next
        print "Time Step " + str(t) + " done."
        
    #File saving
    #print solution[1,-5,-5,:]
    #print solution[tsteps-1,-5,-5,:]
    
    return solution
    
    
def main():
    K_x = spectral(1)
    K_y = spectral(2)
    K_z = spectral(3)
    
    np.save("K_x", K_x)
    np.save("K_y",K_y)
    np.save("K_z",K_z)
    
    
    
    
if __name__ == "__main__":
    main()