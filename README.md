# LatticeCUT

With this project, one may study the effect of the effective electron-electron interaction dervied via a continuous unitary transformation on lattice systems.
The code was used to compute the mean-field data and collective excitation spectra presented in 

Enhanced Superconductivity in Proximity to Peaks in Densities of States
J. Althüser, I. M. Eremin, and G. S. Uhrig
https://doi.org/10.48550/arXiv.2512.11451

and

Secondary Collective Excitations in Intermediate to Strong-Coupling Superconductors
J. Althüser and G. S. Uhrig
https://doi.org/10.48550/arXiv.2605.20059



### Requirements
- C++ 20 and a functioning compiler (tested with g++ 13.3.0 on WSL and g++ 11.5.0 on Red Hat)
- Eigen (tested with version 3.4.1) https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost (tested with version 1.78.0) https://www.boost.org/
- OpenMPI (tested with version 4.1.1) https://www.open-mpi.org/
- OpenMP https://www.openmp.org/
- nlohmann/json.hpp (tested with version 3.11.3) https://github.com/nlohmann/json
- [recommended] CMake 3.30 or newer (tested with version 3.31.8)  https://cmake.org/
- [optional] Intel MKL for BLAS and LAPACK (tested with version 2025.2.1) https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html 

- the mrock library located in `../../PhdUtility/`.
- Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

| Option | Type / Values | Description |
|---|---|---|
| `phonon_coupling` | `<float>` | Effective electron-electron interaction strength *g*. |
| `local_interaction` | `<float>` | Anderson-Morel pseudopotential / Hubbard-like repulsion *U*. |
| `fermi_energy` | `<float>` | Fermi energy. |
| `omega_debye` | `<float>` | Fraction of discretization points $N$ in which the phonon coupling is active. |
| `beta` | `<float>` | Inverse temperature. `-1` is treated as infinity. |
| `dos` | `free_electronsN`, `sc`, `fcc`, `bcc`, `hc`, `single_peakN`, `double_peakN` | Which DOS to use. |
| `N` | `<int>` | Number of abscissae. |
| `output_folder` | `<string>` | Data directory. |
| `investigated_operator` | `0`, `1`, `255` | Operator type for the Green's functions. `0` = full, `1` = near epsilon = 0, `255` = preset defined in the header. |
| `collective_modes` | `<bool>` | Whether to compute the Green's functions. |
| `T_C` | `<bool>` | Whether to compute the critical temperature *T_c*. |
| `complex` | `false` | Whether to use complex numbers, as a sanity check. |



## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine (-march=native).
There are additionally the targets `icelake` and `cascadelake`, which will built for the corresponding CPU architecture.


## Running the program

Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.
such that you have the directory `../commutators/hubbard/` filled with the results of the commutators.
Then, (after building of course), you may create or edit a parameter file in the `params` directory.
Run the program with `./path/to/executable path/to/param/file.config`.
By default, the executable will be located in `./build/default/`.
The result will be saved into `../../data/lattice_cut/<output_folder in the parameter file>/<string generated from the model parameters>/<filename (e.g. resolvents.json.gz)>`.

If full diagonalization (for the operator amplitudes) is request, the build dir will be appended by `_ed`, for example `default_ed` instead of `default`.


## Testing

`make test` will build and run the program with a set of parameters listed in the `tests` directory.
It will then call the plot script in the same dir to visualize the test data.
As before, it requires a previously completed run of `FermionCommute`, which must save the commutators to `../commutators/hubbard/`.
The program will save the simulation data and consequently plot it.

The plots will be saved to `build/default_ed/<plot_name>.pdf`.
For doing so, python with matplotlib, numpy, and pandas is required.
Moreover, the python modules of `../../PhdUtility/python` are needed to evaluate the continued fractions.

### Expected results
#### Standard
Operator amplitudes should exhibit oscillatory behavior with 2 additional roots per mode.

#### Enhanced
Should display subgap peaks (gap is hard to distinguish from the continuum due to the small discretization).
Additionally: Large weight inside the continuum around ω=0.3.