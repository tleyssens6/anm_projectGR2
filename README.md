# anm_projectGR2

This is our smoothed particle hydrodynamics simulation. 
This code simulates the breaking of a dam in 3 different configurations : 
1) simple dam break on a flat and dry surface
2) dam break at a certain altitude into a dry bed
3) dam break at a certain altitude into a pool of fluid beneath

The default simulation is simulation number 3.
To run the different simulations, the user should uncomment the simulation he wishes to see, inside the main function.

The SPH.c file contains most of the "important" functions used in the SPH resolution scheme. 
Particle.c contains the functions linked to the particles and their neighborhood search. 
in kernel.c, the computation of the smoothing function may be found.
in derivatives.c, the different spatial derivatives are computed. 

Credits to : 

Akilan Mathiazhagan
Louis Alsteens
Thomas Leyssens



NB : there seems to be slight unsolvable bug when running the simulation on windows (but not on mac...). If this is the case, please refer to the sph.c function, and comment lines 34 to 38. It will work, but the verlet algorithm will just not be as effective.
