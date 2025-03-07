# A vendor-agnostic representation of MRI pulse sequences

<!-- ![logo](doc/logo.svg) -->

This repository contains MATLAB code for converting a 
[Pulseq](https://pulseq.github.io/) 
MRI pulse sequence file
to a format suitable for GE (and perhaps other) scanners.

For usage information and examples, and for a description
of the Pulseq GE interpreter and its relation to the current repository, 
see: https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2

The idea behind this repository is to represent an MRI pulseq sequence as a 
(typically small) collection of parent/prototype
[Pulseq](https://pulseq.github.io/)
blocks, that are played out repeatedly during the scan with different
RF/gradient amplitudes, RF/DAQ frequency/phase offsets, etc:

![model](doc/model.svg)

The goal of this repository is to provide an **open standard** specification
for encapsulating this sequence description.
At present, this specification is formally expressed as a nested C struct named **Ceq**, however
equivalent representations in other programming languages may be used in practice.


