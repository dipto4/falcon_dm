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

### Usage
The code handles all initial conditions (including global simulation parameters) using a file called ```config.ini```. A sample ```config.ini``` is present in the ```examples/``` directory. This file must be present inside the simulation directory.

A ```config.ini``` looks like the following

```
[Globals]
#0 to use ASCII initial conditions file, 1 to restart from a snapshot_*.hdf5 file
restart = 0

filename = <enter the name of your initial conditions file>
#if using restart = 0, file should be in ASCII format with the columns as M X Y Z VX VY VZ
#IMPORTANT: The first two rows MUST CONTAIN THE IMBH AND THE COMPACT OBJECT RESPECTIVELY

#eta controls the timestep parameter
eta = <enter the timestep parameter here, typical values are 0.025-0.01> 

#soft controls the plummer softening value between the IMBH and the DM particles
soft = <enter the softening length in chosen length units>

#if post newtonian terms are to be used, the following values must be set
use_pn = 1
use_precession = 0
use_radiation = 1

#final time (in G=1 units) and frequency of snapshot (in G=1 units)
t_end = 30.0
dt = 0.01

#if PN terms are used, value of c must be set in units of the base units such that G=1
#here base units are (MSun, AU, day)
c = 10065.3201203821918170

```
The initial conditions file must be either in ASCII format or HDF5 format (for restarts). If restarting, the snapshot should be named as ```snapshot_*.hdf5```.

When using an ASCII file, it is important to note that the data must be in the following format:
```Mass Posx Posy Posz Velx Vely Velz```. Check out the ```examples/``` directory for a few test ICs.

#### Important note:
The first row of the initial conditions file MUST contain the IMBH and the second row MUST contain the compact object. Otherwise the code is going to produce incorrect outputs.

If the code executes properly, the following data should be printed to the standard output every ```dt``` units:
``` CoM_x CoM_y CoM_z TotalEnergy```
The total energy is computed as the total mechanical energy (potential + kinetic) of the system. In case of PN term usage, this is not going to be conserved exactly (using PN precession/radiation terms). This is a result of how the energy is computed rather than any issue with the code/integrator.

The code will output one more line after that which indicates the amount of wall-clock time spent computing the step. 

### Reading the snapshots
The code outputs the snapshots in HDF5 format. The user may choose to read the snapshots in python using ```h5py``` module. For example a sample script may look like the following:

```
import h5py
import numpy as np

data = h5py.File('snapshot_0.hdf5','r')

m = np.array(data['Mass']) # reads from the mass column
x = np.array(data['Posx']) # reads from the Posx data column
y = np.array(data['Posy']) # reads from the Posy data column
z = np.array(data['Posz']) # reads from the Posz data column
vx = np.array(data['Velx']) # reads from the Velx data column
vy = np.array(data['Vely']) # reads from the Vely data column
vz = np.array(data['Velz']) # reads from the Velz data column

# other data attributes are available such as the potential of each particle
pot = np.array(data['Potential'])
```
It is important to note that the order of the particles in the snapshot is preserved. Thus the 3rd particle in the IC file is the same as the 3rd particle in say snapshot 56.

Additionally, if the user wants to read in the simulation time, they should add the following line to their script
```
time = data['Header'].attrs['Time']
```

### Known issues
1. In case of PN term usage, the energy is not going to be conserved exactly (using PN precession/radiation terms). This is a result of how the energy is computed rather than any issue with the code/integrator. We plan on using a different energy calculation criterion to fix this issue.
2. Due to the approximated force calculation, energy/momentum is not going to be conserved exactly. This is again, not an issue with the code/integrator but rather how the energy/momentum is calculated. The code has been validated against other publicly available N-body codes and we find that there are no differences with other codes with regards to the dynamics of the central binary.

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
