__*A molecular quantum dynamics and spectroscopy package*__

# This package is currently capable of the following calculations:

<a href="pldm_info.html"><b>PLDM</b></a> (Partially Linearized Density Matrix)
---
* <a href="calculate__pldm__redmat_8f90.html"><b> Reduced density matrix</b></a> dynamics of a system-bath model where the bath consists of a set of harmonic 
 oscillators.
 
* <a href="calculate__pldm__absorption_8f90.html"><b> Linear optical response function</b></a> for a system-bath model where the bath consists 
of a set of harmonic oscillators.

* <a href="calculate__pldm__nonlinear_8f90.html"><b> Third order optical response function</b></a> for a system-bath model where the bath consists 
of a set of harmonic oscillators. (__currently under construction__)

* <a href="calculate__ipldm__redmat_8f90.html"><b> Iterative reduced density matrix (IPLDM)</b></a> dynamics of a system-bath model where the bath consists of a set of harmonic 
 oscillators. This calculation employs the "focusing" procedure, outlined in the IPLDM example directory.

<a href="twa_info.html"><b>TWA</b></a> (Truncated Wigner Approximation)
---

* <a href="calculate__twa__redmat_8f90.html"><b> Reduced density matrix</b></a> dynamics
of a system-bath model where the bath consists of a set of harmonic oscillators.
(__currently under construction__)

<a href="sqc_info.html"><b>SQC/MM</b></a> (Symmetrical Quasi-Classical / Meyer-Miller)
---
* <a href="calculate__sqc__redmat_8f90.html"><b> Reduced density matrix</b></a> dynamics of a system-bath model where the bath 
consists of a set of harmonic oscillators. Calculation
can be performed with the original *square* windowing scheme or the *triangular* windowing scheme.

* <a href="calculate__sqc__absorption_8f90.html"><b> Linear optical response function</b></a> for a system-bath 
model where the bath consists of a set of harmonic oscillators. Calculation can be performed with the original *square*
windowing scheme or the *triangular* windowing scheme.
(__Currently in submission to the Journal of Chemical Physics__)


<a href="equilibrium_info.html"><b>EQUILIBRIUM REDUCED DENSITY MATRIX</b></a>
---
* <a href="calculate__equilibrium__site_8f90.html"><b> Equilibrium </b></a> Imaginary time path integral calculation of the exact equilibrium 
__reduced density matrix__ in the site basis for a system that has populations bi-linearly
coupled to a bath of harmonic oscillators.

# Compilation

Compilation of this program requires CMake minimum version 3.2. 
To compile this program (from the current directory), execute `mkdir build && cd build` 
followed by `cmake ../ && make`. To run tests, execute `make test`. The executables 
will be in the `run` directory. To build the documentation, execute `make docs`.


__JP gratefully acknowledges the Molecular Sciences Software Institute for funding the development of the MQDS package,
 as well as the mentorship of Dr. Benjamin Pritchard.__
\htmlonly 
<div><a href="http://www.molssi.org"><img src="MolSSI-Logo-2.jpg"/></a></div> 
\endhtmlonly

