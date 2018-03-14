# IPLDM examples

The current examples consist of:
---------------------------------------

* two-state IPLDM reduced density matrix calculation in diabatic basis with harmonic bath and bilinear system-bath coupling
* two-excited-state IPLDM absorption spectrum calculation in diabatic basis with harmonic bath and bilinear system-bath coupling

NOTES
---
* The input file __must__ be processed using `python process_input.py` before running a calculation.
* This calculation employs the "focusing" procedure which consists of an importance sampling of the density
matrix elements as well as a "steepest descent" analysis of the intermediate integrals over mapping 
variable initial conditions.