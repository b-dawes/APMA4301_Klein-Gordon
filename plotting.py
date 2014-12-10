# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:45:11 2014

@author: tfuser
"""
import numpy as np
import matplotlib.pyplot as plt
'''
# Takes in k-matrix and solves for value at t,x,y,z
def F(K,k_max,dk,t,x,y,z):
    out = 0
    for k_x in xrange(-k_max,k_max,dk):
        for k_y in xrange(-k_max,k_max,dk):
            for k_z in xrange(-k_max,k_max,dk):
                out = out + K[t][k_x][k_y][k_z]*np.exp(1j*(k_x*x+k_y*y+k_z*z))
    return out

# Takes in k-matrix and converts to position space matrix
def toPosition(K,k_max,dk,x,y,z):
    X = np.zeros((K.shape[0],len(x),len(y),len(z)))
    for t in xrange(0,1):
        for i in xrange(0,len(x)):
            for j in xrange(0,len(y)):
                for k in xrange(0,len(z)):
                    print 't=',t,' i=',i,' j=',j,' k=',k
                    X[t][i][j][k]=F(K,k_max,dk,t,x[i],y[j],z[k])
    return X
'''

def toPosition(K):
    np.fft.ifftn(K)
# Takes in X matrix and plots at time t on plane where the dimension d = val
def plotPlane(X,t,d,val):
    plt.figure()    
    if d=='x':
        plt.pcolor(X[t][val][:][:])
        plt.xlabel('y')
        plt.ylabel('z')
    if d=='y':
        plt.pcolor(X[t][:][val][:])
        plt.xlabel('x')
        plt.ylabel('z')
    if d=='z':
        plt.pcolor(X[t][:][:][val])
        plt.xlabel('x')
        plt.ylabel('y')
    plt.colorbar()
    plt.title('')
    plt.draw()
    plt.show()
    print 1
                        
def main():
    print 'Importing variables...'
    K_x = np.load('K_x.npy')
    #K_y = np.load('K_y.npy')
    #K_z = np.load('K_z.npy')
    
    params = np.load('params.npy')
    
    print 'Converting to position space...'
    #X_x = toPosition(K_x,k_max,dk,x,x,x)
    print K_x[0][5][5][:]
    X_x = np.zeros(K_x.shape)    
    for i in xrange(0,K_x.shape[0]-1):   
        X_temp = np.fft.ifftn(K_x[i][:][:][:])
        X_x[i][:][:][:] = X_temp
    
    #print X_x[:][:][6]
    print 'Plotting...'
    for i in xrange(X_x.shape[0]-1):
        plotPlane(X_x,X_x.shape[0]-1-i,'z',6)
    
if __name__ == "__main__":
    main()