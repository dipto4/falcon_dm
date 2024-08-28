## Falcon-DM

### Description
Lightweight N-body code written in C++, specifically tuned for simulations of Intermediate Mass Ratio Inspirals (IMRIs) embedded in dark matter (DM) spikes.

### Related publications
[_Examining the Effects of Dark Matter Spikes on Eccentric Intermediate Mass Ratio Inspirals Using  N-body Simulations_(2024) **Mukherjee, D.**, Holgado, A. Miguel, Ogiya, Go & Trac, H. Monthly Notices of the Royal Astronomial Society, 533(2) 2335-2355](https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/stae1989/7737663)

### Features
1. 2nd order Drift-Kick-Drift integrator using the symplectic HOLD scheme
2. Symmetrized, individual, time-steps for accurate time-integration
3. Post-Newtonian (PN) effects upto PN2.5 using the auxiliary velocity algorithm

### Additional planned features
1. Inclusion of a 4th order symplectic integrator
2. Inclusion of PN terms upto 3.5 order
3. Enhanced vectorization using SIMD 

## Instructions

### Requirements
1. Compiler supporting C++17 standards
2. Boost library (for reading .ini files)
3. HDF5 library (for reading and writing snapshots)

### Installation
1. Navigate to the directory and modify the Makefile according to the compiler suite present on your workstation along with any compiler flags you want. Note that ```-fopenmp``` flags must be enabled as the ```OpenMP``` library is used to calculate the time spent during integration (other than the energy calculation parallelization
2. Ensure that the relevant directories for the HDF5 and Boost libraries are included properly in the Makefile
3. Type ```make```. This should create an executable called ```falcon```. 
4. The code handles all initial conditions (including global simulation parameters) using a file called ```config.ini```. A sample ```config.ini``` is present in the ```examples/``` directory. This file must be present inside the simulation directory.

If you use this code, please cite Mukherjee et al. (2024) 

```
@ARTICLE{2024MNRASMukherjee_FalconDM,
       author = {{Mukherjee}, Diptajyoti and {Holgado}, A. Miguel and {Ogiya}, Go and {Trac}, Hy},
        title = "{Examining the effects of Dark Matter spikes on eccentric Intermediate Mass Ratio Inspirals using N-body simulations}",
      journal = {\mnras},
         year = 2024,
        month = aug,
          doi = {10.1093/mnras/stae1989},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.tmp.1946M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```
