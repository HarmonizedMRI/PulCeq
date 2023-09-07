function seq2ge(seqarg, ofname, varargin)
% Convert a Pulseq file (http://pulseq.github.io/) to a .tar file for GE scanners
%
% Inputs
%   seqarg = a seq object, or name of a .seq file
%
% Ouputs
%   none
%
% Effect
%   Writes a .tar file containing all the required scan files to run on GE scanners

% Parse inputs
% Read in system specs defined from .seq file
seq = mr.Sequence(); seq.read(seqarg); sys = seq.sys;

% Define sysGE
arg.sysGE = toppe.systemspecs(...
    'psd_rf_wait', 200,... % us
    'psd_grd_wait', 200,... % us
    'maxGrad', sys.maxGrad/sys.gamma*100,... % G/cm
    'maxSlew', 1.01*sys.maxSlew/sys.gamma/10,... % G/cm/ms. Factor > 1 is fudge factor to avoid exceeding limit after interpolating to GE raster time.
    'rfDeadTime', sys.rfDeadTime*1e6,... % us
    'rfRingdownTime', sys.rfRingdownTime*1e6,... % us
    'adcDeadTime', sys.adcDeadTime*1e6...  % us
);

% Use explicitly provided sysGE (function from MIRT toolbox)
arg = vararg_pair(arg, varargin);

% Convert the .seq file/object to the PulCeq representation
ceq = seq2ceq(seqarg);

% Write to TOPPE files
ceq2ge(ceq, arg.sysGE, ofname);
