
![logo](doc/logo.svg)

**[NB! We're just getting started -- current development is happening in the 'dev' branch]**

A C representation of arbitrary magnetic resonance pulse sequences, 
based on a "scaled parent block" description.

The idea is to represent the sequence as a (typically small) collection of parent/prototype
[Pulseq](https://pulseq.github.io/)
blocks, that are played out repeatedly during the scan with different
RF/gradient amplitudes, RF/DAQ frequency/phase offsets, etc:

![model](doc/model.svg)

The goal of this repository is to provide an **open standard** specification
for encapsulating this sequence description in C code.

This specification is contained in the ./src/ folder, which contains the following:
* **pulCeq.h**: defines `Ceq`, a nested struct containing the entire sequence, and
function declarations (interfaces) for various tasks including
   * serializing to file
   * allocating and freeing up memory
   * <...>
* **pulCeq.c**: implementation

## Example usage

* Planned: PulCeq.h will be used in the upcoming version of the Pulseq sequence interpreter for GE scanners (TOPPE v6)

