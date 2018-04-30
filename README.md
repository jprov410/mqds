[![Build Status](https://travis-ci.org/jprov410/mqds.svg?branch=master)](https://travis-ci.org/jprov410/mqds)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![codecov](https://codecov.io/gh/jprov410/mqds/branch/master/graph/badge.svg)](https://codecov.io/gh/jprov410/mqds)

![alt text](https://github.com/jprov410/mqds/blob/master/images/logo.jpg)

__*A molecular quantum dynamics and spectroscopy package*__

#### This package is currently capable of the following calculations:


* __PLDM__ (Partially Linearized Density Matrix) calculation of the __reduced density
 matrix__ of a system-bath model where the bath consists of a set of harmonic 
 oscillators.

* __PLDM__ (Partially Linearized Density Matrix) calculation of the __linear 
optical response function__ for a system-bath model where the bath consists 
of a set of harmonic oscillators.

* __PLDM__ (Partially Linearized Density Matrix) calculation of the __third order 
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

* __SQC__ (Symmetrical Quasi-Classical) calculation of the __linear optical response function__ 
of a system-bath model where the bath consists of a set of harmonic oscillators.

* __EQUILIBRIUM__ Imaginary time path integral calculation of the exact equilibrium 
__reduced density matrix__ in the site basis for a system that has populations bi-linearly
coupled to a bath of harmonic oscillators.

More Information
----

* __Travis CI build status:__ gfortran - 4.9, 5, 6, 7 & pgf90
[![Build Status](https://travis-ci.org/jprov410/mqds.svg?branch=master)](https://travis-ci.org/jprov410/mqds)
* __Codecov report:__ coverage reported by gfortran-6 build on Travis 
[![codecov](https://codecov.io/gh/jprov410/mqds/branch/master/graph/badge.svg)](https://codecov.io/gh/jprov410/mqds)

* __Documentation:__ Documentation can be found at ( https://jprov410.github.io/mqds ) 

## Compilation

Compilation of this program requires CMake minimum version 3.2. 
To compile this program (from the current directory), execute `mkdir build && cd build` 
followed by `cmake ../ && make`. To run tests, execute `make test`. The executables 
will be in the `run` directory. To build the documentation, execute `make docs`.

__JP gratefully acknowledges the Molecular Sciences Software Institute for funding the development of the MQDS package__

[![alt text](https://github.com/jprov410/mqds/blob/master/images/MolSSI-Logo-2.jpg)](http://molssi.org)