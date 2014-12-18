This directory contains the files for Brian Dawes and Stephen Carr's project for APMA43001 on the Klein Gordon equation for a massive photon.

The python script kgsimulate.py is used to solve the equation in k-space using a spectral method. The result is output in the .npy format.
The python script plotting.py takes in these .npy file and converts and plots the solution in the xy-plane at z=0 in position-space. It also produces animated gifs of the electric field components and the energy distribution.

Each simulation is contained within a separate subfolder. The subfolder contains the scripts used to run the simulation and the resulting gifs.

In the file name, mX means that the simulation ran with m=X. sigmoid indicates that sigmoid dampening was applied in the simulation.

steadystate_tests contains tests involving different dampening filters. The file conv_plots.py is used to plot the residual of the different filters.

paper contains the tex source, images, and the output pdf for our report.