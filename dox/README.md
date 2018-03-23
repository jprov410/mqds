__*A molecular quantum dynamics and spectroscopy package*__

#### This package is currently capable of the following calculations:


* __PLDM__ (Partially Linearized Density Matrix) calculation of the __reduced density
 matrix__ of a system-bath model where the bath consists of a set of harmonic 
 oscillators.
 


* __PLDM__ (Partially Linearized Density Matrix) calculation of the __linear 
optical response function__ for a system-bath model where the bath consists 
of a set of harmonic oscillators.

* __IPLDM__ (Iterative Partially Linearized Density Matrix) calculation of the __reduced density
 matrix__ of a system-bath model where the bath consists of a set of harmonic 
 oscillators. This calculation employs the "focusing" procedure, outlined in the IPLDM example directory.

* __TWA__ (Truncated Wigner Approximation) calculation of the __reduced density matrix__ 
of a system-bath model where the bath consists of a set of harmonic oscillators.
(__currently under construction__)

* __SQC__ (Symmetrical Quasi-Classical) calculation of the __reduced density matrix__ 
of a system-bath model where the bath consists of a set of harmonic oscillators. 

## Compilation

Compilation of this program requires CMake minimum version 3.2. 
To compile this program (from the current directory), execute `mkdir build && cd build` 
followed by `cmake ../ && make`. To run tests, execute `make test`. The executables 
will be in the `run` directory. To build the documentation, execute `make docs`.