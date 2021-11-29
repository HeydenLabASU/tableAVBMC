# tableAVBMC
Aggregation-Volume-Bias Monte Carlo simulation with tabulated potentials

This single file program can be compiled with any ANSI C compiler, 
e.g. on most UNIX platforms via:

gcc tableAVBMC.c -o tableAVBMC
  
The purpose of this program is to run an agregation-volume-bias Monte Carlo (AVBMC) 
simulation in the canonical ensemble with interacting classical particles in 
3 dimensions with periodic boundary conditions. 

Simulation parameters are provided in a single input file, which among other 
simulation parameters specifies one or multiple files with tabulated potentials 
for pair-wise additive interactions between particles.

Up to 10 distinct particle types are supported and require tabulated potentials 
for each possible interaction pair, i.e. for n distinct particle types you'll 
need (1/2)*n*(n+1) potentials.

The AVBMC algorithm requires defining an in and an out region, which is implemented 
by a particle-independent minimum and maximum distance between particle centers.

An exlucsion radius can be defined, which generates an infinite energy if any 2 
particles are found within a distance smaller than the exclusion radius.

A commented input file (template.input) is provided, which provides a usage example 
for a simulation with a single particle type (explanations are given on how to use 
multiple particle types). Comments in the output file (lines starting with "#") 
are automatically discarded by the program and do not need to be removed to run 
the program.

If the program (tableAVBMC) and the template input file (template.input) are 
located in the same directory, simple run the program on the command line with:

./tableAVBMC template.input

For details on the AVBMC algortithm(s), pleasr read and cite these articles:
References:
Bin Chen and J. Ilja Siepmann
A Novel Monte Carlo Algorithm for Simulating Strongly Associating Fluids:  
Applications to Water, Hydrogen Fluoride, and Acetic Acid
J. Phys. Chem. B 2000, 104, 36, 8725–8734
https://doi.org/10.1021/jp001952u

Bin Chen and J. Ilja Siepmann
Improving the Efficiency of the Aggregation−Volume−Bias Monte Carlo Algorithm
J. Phys. Chem. B 2001, 105, 45, 11275–11282
https://doi.org/10.1021/jp012209k

