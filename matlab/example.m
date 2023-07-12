% create a .seq file
addpath ~/github/HarmonizedMRI/B0shimming/sequence/Pulseq/
writeB0;
fn = 'b0';

% convert to PulCeq representation
addpath ~/github/HarmonizedMRI/PulCeq/matlab/ 
fprintf('Converting to PulCeq struct:\n')
ceq = seq2ceq([fn '.seq'], 'verbose', false);

% Write to TOPPE files
addpath ~/github/toppeMRI/toppe/
sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'maxSlew', 20, ...                % G/cm/ms
    'maxRF', 0.15, ...                % Gauss. Must be >= peak RF in sequence.
    'maxSlice', 60, ...
    'rfDeadTime', 72, ...             % us
    'rfRingdownTime', 54, ...         % us
    'psd_rf_wait', 148, ...           % UHP @ RX28 = 148us
    'psd_grd_wait', 156, ...          % UHP @ RX28 = 156us
    'segmentRingdownTime', 104, ...   % us
    'adcDeadTime', 40);               % us

fprintf('Writing TOPPE files:\n');
ceq2ge(ceq, sysGE, [fn '.tar'], true);

% The file 'mprage.tar' can now be executed on a GE scanner
