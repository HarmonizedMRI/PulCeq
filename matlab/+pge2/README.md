# Functions for checking and validating a 'Ceq' sequence representation

This folder contains the `+pge2` namespace, which defines functions for
checking and plotting a Ceq object before writing it to a  .pge file with `writeceq.m`,
and **validating** the GE interpreter output (via WTools or scanner) against the original .seq file.

The latter is based on GE's 
[WTools simulator](https://github.com/jfnielsen/TOPPEpsdSourceCode/blob/UserGuide/v7/simulate.md) 
and the accompanying 'Pulse View' plotter
which is the best way to test and visualize a .pge file before going to the scanner.


Workflow example:
```matlab
% Load the .seq file and convert to a Ceq object
fn = 'gre2d';        % .seq file name
ceq = seq2ceq([fn '.seq']);

% Check the ceq object:
% First define hardware parameters.
psd_rf_wait = 100e-6;  % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6; % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
coil = 'xrm';          % MR750. See pge2.opts()
sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

% Check if 'ceq' is compatible with the parameters in 'sysGE'.
% Specifically, this step checks PNS and RF/gradient limits.
pars = pge2.check(ceq, sysGE);

% Plot part of the sequence 
S = pge2.plot(ceq, sysGE, 'timeRange', [0 0.02], 'rotate', false); 

% Write ceq object to file
% pislquant is the number of ADC events used to set Rx gains in Auto Prescan
writeceq(ceq, [ fn '.pge'], 'pislquant', 10);

% After simulating in WTools/VM or scanning, grab the xml files
% and compare with the seq object, e.g.:
warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals
xmlPath = '~/transfer/xml/';
seq = mr.Sequence();
seq.read([fn '.seq']);
pge2.validate(ceq, sysGE, seq, xmlPath, 'row', 'all', 'plot', true);
```

Also see https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2
for further details and the most up to date information.

