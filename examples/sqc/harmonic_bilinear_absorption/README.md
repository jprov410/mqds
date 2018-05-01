# PLDM reduced density matrix 2-state, harmonic bath, bilinear coupling

* __run.in:__ Input file where user can define specific values for input variables
* __keyword_library.py:__ Module that contains all keywords as well as their default values.
* __hel.in:__ Contains the matrix elements for the system part of the hamiltonian in the diabatic basis
* __continuumsd.in:__ Contains the spectral density that defines the frequency-dependent system-bath coupling strength. Format: First column is frequency, remaining are J(w).
* __couplings.in:__ Determines the number of baths and which states are coupled to each bath.
* __dipole.in:__ Contains the matrix elements of the dipole operator

NOTES
---
* Be sure to include the number of states `nstates` in the input file to reflect the fact that there is now a ground electronic state included in the calculation.
* The input hamiltonian __must__ have elements for the ground electronic state which is not coupled to the excited state manifold.
* The `dipole.in` file should only have nonzero elements that correspond to transitions between the ground and excited state manifold (assuming that there is no permanent dipole).

METHODOLOGY DETAILED IN:
---
* __IN SUBMISSION. 2018__ 