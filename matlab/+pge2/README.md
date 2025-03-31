# Functions for checking and inspecting a 'Ceq' sequence representation

## pge2.validate()

Check a .pge file against GE MRI scanner hardware requirements.
Example:

```matlab
% Convert a Pulseq file to a Ceq sequence object (a nested struct)
ceq = seq2ceq('gre2d.seq');

% Define hardware parameters
psd_rf_wait = 150e-6;  % RF events are delayed by this amount to account for gradient delay (s)
psd_grd_wait = 120e-6; % ADC window is delayed by this amount to account for gradient delay (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
sys = pge2.getsys(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, gamma);

% Check if 'ceq' is compatible with the parameters in 'sys'
pge2.validate(ceq, sys);

% Modify ceq until the above call returns no errors

% Write sequence to file
writeceq(ceq, 'gre2d.pge');
```

As usual, see 
https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2
for further details and the most up to date information.

## Plot a base segment

To plot the first base segment from the above example, do:
```matlab
S = pge2.constructvirtualsegment(ceq.segments(1).blockIDs, ceq.parentBlocks, sys, true);
```
Note that all waveforms are shown with normalized (positive) amplitude.


## Units

Units used here is a mix of Pulseq and GE conventions, and are as follows:
```
RF:                  Gauss (typical GE hardware limit is 0.25 G)
Gradients:           Gauss/cm    
Slew rate:           Gauss/cm/ms
Times:               sec
gyromagnetic ratio:  Hz/Gauss  (4257.6 Hz/Gauss for 1H)
```
