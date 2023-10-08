# Interactive topology optimization notebook with overhang filter
by Emiel van de Ven

This repository is supplementary material to "Overhang control in topology optimization: a comparison of continuous front propagation-based and discrete layer-by-layer overhang control", E. van de Ven et al, 2020. It contains a jupyter notebook that can be run by starting jupyter notebook with:
```
jupyter notebook
```
and opening `TopOptAMapplet.ipynb`.

The topolgoy optimization is based on the python implementation of the [88-line topolgoy optimization code by Niels Aage and Villads Egede Johansen](http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python), with some minor modifications to incorperate the overhang filters. The layer-by-layer filter that is implemented is a python adaption of the supplementary material of ["An additive manufacturing filter for topology optimization of print-ready designs", M. Langelaar (2016)](https://link.springer.com/article/10.1007/s00158-016-1522-2). The front-propagation based overhang filter is based on ["Continuous front propagation-based overhang control for topology optimization with additive manufacturing", E van de Ven et al.](https://link.springer.com/article/10.1007/s00158-017-1880-4), while the improved version is described in the current publication.

The front propagation code is written in c++: because of the sequential nature of front-propagation, the code cannot be vectorized easily, and is better written in a compiled language. The repository contains compiled versions for windows and unix systems. If your system is not supported, recompile FPAMFilter.cpp.

For windows:
```
cl /LD /Ox /Qpar FPAMFilter.cpp
```
For Unix:
```
g++ -shared -fPIC -O3 -o libFPAMFilter.so FPAMFilter.cpp
```

The file `FPAMFilter.cpp` defines a class `FPAMFilter`, which has four functions: 
* `FPAMFilter`
* `evaluate`
* `sens`
* `updateArrivalTime`

The constructor function `FPAMFilter` does some initialization of variables. The `evaluate` function takes a density field and returns the overhang filtered density field. `sens` performs the adjoinst sensitivity calculation. `updateArrivalTime` calculates the arrival time of an element given two adjacent elements with known arrival times.

The `evaluate` functions is structured as follows:
* Initialize the arrival times of the nodes on the base-plate, set status of bottom nodes to Accepted, rest to Candidate.
* Start the front propagation loop. In each loop:
  1. Find the Candidate element `emin` with the lowest arrival time, set status to Accepted.
  2. Loop over the Candidate neighbours of `emin`, update their arrival times.
  3. Repeat until all nodes are Accepeted.
* Calculate the densities from the arrival time field by calculating the delay field.

The `sens` function simply loops over all the nodes in opposite order as in which they are Accepeted in the `evaluate` function, and calculates the adjoint sensitivites.

In `updateArrivalTime`, given the arrival times at two locations, the arrival time at a third location is calculated according to the speed function chosen.

The file `FPAMFilter.cpp` is well commented, please check it out for further details.