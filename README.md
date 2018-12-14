#### Available data and codes for the work: 
"Chromatic" Neuronal Jamming in a Primitive Brain

### Contents in each folder are described below.

#### exp_data/
experimental data of 34 planarian samples of 6 imaged neuropeptides (in name): matrix of neuronal cell centroids 
experimental data of 8 double-FISH planarian samples of 3 neuropeptides  

#### samples_lammps_run/
sample simulation data 
each folder contains results of sample LAMMPS runs. Folder name indicate controllinig parameters of color number, packing density, and homotypic interaction cut-off length. 

nC: number of unique colors in the system  

rho: packing fraction

rC: homotypic cut-off distance in terms of heterotypic interaction range

In each folder, 

in.#.min are lammps input files

log.#.out are lammps output files. Final energy after the minimization can be found in it.

min.#.atom are minized packing structures

cluster.#.dat are clustering results. Each line stores cluster information for beads of a single color. For example, the 1st line is clustering information for beads of type 1 in the configuration. The sequence of the beads are the same as they appear in min.#.atom

#### analysis_codes/
contains code for analysis of experimental data and simulation data

experiments/: MATLAB analysis of centroids of neuronal cells including Voronoi tesselation and fitting using Granocentric model   

simulation/: analysis of clustering

  The cluster.py will read in LAMMPS output (for example min.1.atom in sample run folders) and conduct clustering for beads of each color. To construct the distance matrix distance, it will use a subroutine based on c (dist.c) to speed up. To use it, you have to recompile dist.c and change the path to the location of "dist" subroutine command in cluster.py. You may also need to change the location where you store the LAMMPS output files.




