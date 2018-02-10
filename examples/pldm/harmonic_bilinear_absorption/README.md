# PLDM reduced density matrix 2-state, harmonic bath, bilinear coupling

* __run.in:__ Raw input file where user can define specific values for input variables
* __process_input.py:__ Python script that produces the input file that the __mqds__ executable interprets. This file contains values for _all_ input variables. If values are not provided in the raw __run.in__ file, they are assigned their default values. Execute `python process_input.py` to produce the __processed_run.in__ input file.
* __keyword_library.py:__ Module that contains all keywords as well as their default values.
* __hel.in:__ Contains the matrix elements for the system part of the hamiltonian in the diabatic basis
* __continuumsd.in:__ Contains the spectral density that defines the frequency-dependent system-bath coupling strength. Format: First column is frequency, remaining are J(w).
* __couplings.in:__ Determines the number of baths and which states are coupled to each bath.