% create a .seq file
addpath ~/github/HarmonizedMRI/MP-RAGE/sequence/Pulseq/
% writeMPRAGE;

% convert to PulCeq representation
addpath ~/github/HarmonizedMRI/PulCeq/matlab/ 
%ceq = seq2ceq('mprage.seq');

% Write to TOPPE files
addpath ~/github/toppeMRI/toppe/
sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'B0', 3.0, ...
    'maxSlew', 20, ...                % G/cm/ms
    'maxRF', 0.15, ...                % Gauss. Must be >= peak RF in sequence.
    'timessi', 100, ...               % us
    'rfDeadTime', 72, ...             % us
    'rfRingdownTime', 54, ...         % us
    'segmentRingdownTime', 104, ...   % us
    'adcDeadTime', 40);               % us

ceq2ge(ceq, sysGE, 'mprage.tar', true);

% The file 'mprage.tar' can now be executed on a GE scanner
