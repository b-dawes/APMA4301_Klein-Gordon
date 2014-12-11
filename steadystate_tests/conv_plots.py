# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:45:11 2014

@author: tfuser
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
    
def main():
    print 'Importing variables...'
    undamp = np.load('undampened_res_conv.npy')
    sigmoid = np.load('sigmoid_res_conv.npy')
    linear = np.load('linear_res_conv.npy')
    undamp_shift = np.zeros(linear.shape)    
    
    for i in xrange(0,undamp.shape[0]-15):
        undamp_shift[i+15] = undamp[i]
        
    plt.hold(True)
    plt.semilogy(undamp_shift,'r',label='Undampened')
    plt.semilogy(linear,'g',label='Linear')
    plt.semilogy(sigmoid,'b',label='Sigmoid')
    plt.xlabel('time step t')
    plt.ylabel('L2 residual error')
    plt.title('residual convergence for different dampening schemes')
    plt.legend()
        
            
if __name__ == "__main__":
    main()