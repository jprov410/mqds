# TWA examples

The current examples consist of:
---------------------------------------

* two-state sqc reduced density matrix calculation in diabatic basis with harmonic bath and bilinear system-bath coupling

NOTES
---
* The input file __must__ be processed using `python process_input.py` before running a calculation.
* For SQC dynamics, two parameters must be defined in the input file. The variable __zpe__ defines the zero point energy of mapping oscillators (if not defined will take default value of 0.5). __window__ defines the width of rectangular windows and the default value is also 0.5.
* __zpe__ and __window__ __MUST__ maintain positive values.