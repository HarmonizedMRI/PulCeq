Functions for simulating and checking a 'Ceq' sequence representation on GE hardware.

The goal here is to realize the following workflow:
```matlab
% Convert a Pulseq file to a Ceq sequence object (a nested struct)
ceq = seq2ceq('gre2d.seq');

% Define hardware parameters
psd_rf_wait = 150e-6;  % RF pulses are delayed by this amount to account for gradient delay (s)
psd_grd_wait = 120e-6; % ADC window is delayed by this amount to account for gradient delay (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
sys = getsys(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, gamma);

% Check if 'ceq' is compatible with the parameters in 'sys'
pge2.validate(ceq, sys);

% Write sequence to file
writeceq(ceq, 'gre2d.pge');
```

Units used here is a mix of Pulseq and GE conventions, and are as follows:
```
RF:                  Gauss (typical GE hardware limit is 0.25 G)
Gradients:           Gauss/cm    
Slew rate:           Gauss/cm/ms
Times:               sec
gyromagnetic ratio:  Hz/Gauss  (4257.6 Hz/Gauss for 1H)
```
