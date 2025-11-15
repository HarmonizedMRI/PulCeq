function writeentryfile(n, seqname, varargin)
% function writeentryfile(n, seqname, varargin)
% 
% Inputs
%   n        [1]         entry file number (CV1)
%   seqname  string      .seq file name (with or without .seq/.pge extension)
%
% Input options
%   path     string      .pge file location on scanner/WTools simulator

% defaults
arg.path = '/srv/nfs/psd/usr/psd/pulseq/v7/sequences/';

% substitute with provided keyword arguments
arg = vararg_pair(arg, varargin);   % in ../

% strip .seq or .pge extension if present
seqname = replace(seqname, {'.seq', '.pge'}, '');

% write .entry file
fid = fopen([arg.path 'pge' num2str(n) '.entry'], 'wt');

fprintf(fid, '1\n');
fprintf(fid, '%s\n', [arg.path seqname '.pge']);

fclose(fid);
