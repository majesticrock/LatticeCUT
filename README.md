# ContinuumSystem

Code used to compute the mean-field data and collective excitation spectra presented in 

TBD



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

### Effective electron-electron interaction strength *g*
`phonon_coupling <float>`
### Anderson-Morel pseudopotenial *mu-star*
`local_interaction <float>`
### Fermi energy
`fermi_energy <float>`
### Fraction of discretization points (N) in which the phonon_coupling is active
`omega_debye <float>`
### Inverse temperature, -1 is treated as infinity
`beta <float>`
### Which DOS?
`dos <free_electronsN|sc|fcc|bcc|hc|single_peakN|double_peakN>`
### Number of abscissae
`N <int>`
### Data dir
`output_folder <string>`
### Which type of operator to use for the Green's functions. 0 = full, 1 = near epsilon=0
`investigated_operator <0|1>`
### Whether to compute the Green's functions
`collective_modes <bool>`
### Whether to compute the critical temperature *T_c*
`T_C <bool>`

### Required externals
- Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page or https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost https://www.boost.org/
- (optional) Intel MKL https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html - if present on the system, it will be used for Eigen. If not, the standard Eigen routines are employed.
- OpenMPI https://www.open-mpi.org/
- CMake https://cmake.org/
- nlohmann/json.hpp https://github.com/nlohmann/json

### Required internals

- mrock/utility
- mrock/symbolic_operators

see https://github.com/majesticrock/PhdUtility


## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine.
Specifying `cascadelake` or `icelake` will built for the specific CPU architecture.
These two are present on the compute cluster used to generate the data.


## Running the program

Before executing the program, make sure to build and run FermionCommute
https://github.com/majesticrock/FermionCommute/
such that you have the directory `../commutators/lattice_cut/` filled with the results of the commutators.
`./build/latticecut <parameter file>` will run the program using the specified parameter file.
For large scale computations, SLURM scripts are provided.
A few additional bash scripts are provided for ease of running multiple jobs.