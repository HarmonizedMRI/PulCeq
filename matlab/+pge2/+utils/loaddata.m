function dat = loaddata(fn, nADC)
% loaddata - Load GE raw data (ScanArchive) file
%
% function d = loaddata(fn, nADCevents)
%
% Inputs
%   fn     string     ScanArchive file name
%   nADC   [1]        Number of ADC events to load
% Output
%   d      cell array of length nADC (complex)
%
% Usage example:
%    >> addpath ~/Programs/orchestra-sdk-2.1-1.matlab/
%    >> d = pge2.utils.loaddata('mydata.h5', 256);

% Open file
archive = GERecon('Archive.Load', fn);

% get total number of shots
nShots = archive.FrameCount;  
if nargin < 2
    nADC = nShots;
end
assert(nShots >= nADC, sprintf('Data file contains only %d ADC events', nShots));

% read data
clear dat
fprintf(sprintf('Loading %s\n', fn));
textprogressbar('');
for n = 1:nADC
    currentControl = GERecon('Archive.Next', archive);
    dat{n} = currentControl.Data;
    s(n) = max(real(currentControl.Data(:)));
    textprogressbar(n/nADC*100);
end
textprogressbar('');

fprintf('max(real(signal)) = %d\n', max(s));

% save dat dat -v7.3

% [max(real(dat.echo1(:))) max(real(dat.echo2(:)))]
