# Functions for checking and validating a 'Ceq' sequence representation

This folder contains the `+pge2` namespace, which defines functions for
checking and plotting a Ceq object before writing it to a  .pge file with `writeceq.m`,
and **validating** the GE interpreter output (via WTools or scanner) against the original .seq file.

The latter is based on GE's 
[WTools simulator](https://github.com/jfnielsen/TOPPEpsdSourceCode/blob/UserGuide/v7/simulate.md) 
and the accompanying 'Pulse View' plotter
which is the best way to test and visualize a .pge file before going to the scanner.

For the most up to date example of how to use these functions, 
see our 'official' [2DGRE sequence example](https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2/2DGRE).

See https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2
for further details.

