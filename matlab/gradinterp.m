function wavOut = gradinterp(g, rasterIn, rasterOut)
% function wavOut = gradinterp(g, rasterIn, rasterOut)
%
% Interpolate arbitrary gradient to uniform raster
%
% Intputs:
%  g             Pulseq (arbitrary) gradient struct
%                g.waveform     gradient waveform samples (Hz/m)
%                g.tt           sample times (sec)
%                g.first        waveform value at starting edge of 1st raster time
%                g.last         waveform value at ending edge of last raster time
%  rasterIn      input gradient raster time (sec)
%  rasterOut     output gradient raster time (sec)
%
% Output:
%  gout          gradient on uniform raster (Hz/m)

if ischar(g)
    sub_test();
    return;
end

% Round sample times to nearest 100ns to avoid numerical precision error in mod() below
dt = 0.1e-6;
g.tt = round(g.tt/dt)*dt;
rasterIn = round(rasterIn/dt)*dt;
rasterOut = round(rasterOut/dt)*dt;

% initial values
ttIn = g.tt(:);
wav = g.waveform(:);

% If first sample is not at t=0, add g.first.
% This is probably due to g.tt(1) = rasterIn/2
if g.tt(1) > 0
    ttIn = [0; ttIn];
    wav = [g.first; wav];
end

% If last sample is not at end edge of last raster point, add g.last.
if mod(g.tt(end), rasterIn)
    ttIn = [ttIn; ttIn(end) + rasterIn - mod(ttIn(end), rasterIn)];
    wav = [wav; g.last];
end

% Slight distortion that we may have to ask for forgiveness for later: 
% stretch time so end sample is on rasterOut boundary
timeStretchFactor = (ceil(ttIn(end)/rasterOut)*rasterOut) / ttIn(end);
ttIn = ttIn * timeStretchFactor;
dt = 1e-9;
ttIn = round(ttIn/dt)*dt; 
ttOut = [0:rasterOut:ttIn(end)];

wavOut = interp1(ttIn, wav, ttOut); %, 'interp', 'extrap');

return

function sub_test

rasterIn = 10e-6;
g.waveform = [1 1 0]*1e3;       
g.tt = [0 100 400]*1e-6;
g.first = g.waveform(1);
g.last = g.waveform(end);

rasterOut = 4e-6;

gout = gradinterp(g, rasterIn, rasterOut);

tt = (1:length(gout))*rasterOut;

plot(tt, gout, 'bo-');


