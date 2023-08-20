function wavOut = gradinterp(g, rasterIn, rasterOut)
% function wavOut = gradinterp(g, rasterIn, rasterOut)
%
% Interpolate arbitrary gradient to uniform raster
%
% Intputs:
%  g             Pulseq (arbitrary) gradient struct. Must contain the following fields:
%                g.waveform     gradient waveform samples (a.u.)
%                g.tt           sample times (sec)
%                g.first        waveform value at starting edge of 1st raster time
%                g.last         waveform value at ending edge of last raster time
%  rasterIn      input gradient raster time (sec)
%  rasterOut     output gradient raster time (sec)
%
% Output:
%  gout          gradient on uniform raster
%
% To run test function:
%  >> gradinterp('test');

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

% If first sample is not at t=0, add g.first
% This is probably due to g.tt(1) = rasterIn/2
if g.tt(1) > 0
    ttIn = [0; ttIn];
    wav = [g.first; wav];
end

% If last sample is not at end edge of last raster point, add g.last
if mod(g.tt(end), rasterIn)
    ttIn = [ttIn; ttIn(end) + rasterIn - mod(ttIn(end), rasterIn)];
    wav = [wav; g.last];
end

% Output sample times on regular raster.
dt = 1e-9;
ttIn = round(ttIn/dt)*dt; 
ttOutLast = ceil(ttIn(end)/rasterOut)*rasterOut; % add sample at end if necessary to avoid extrapolation
ttOut = [0:rasterOut:ttOutLast];

wavOut = interp1(ttIn, wav, ttOut); %, 'interp', 'extrap');

return


function sub_test

% Define test waveform and sample times.
% Sample times at center of rasterIn intervals, and not on 4us boundary
rasterIn = 10e-6;   % Siemens gradient raster
g.waveform = [1 1 0];   % gradient waveform, a.u.
g.tt = [5 105 205]*1e-6;    % sample times, sec
g.first = g.waveform(1);    % waveform at start edge of first raster time
g.last = g.waveform(end);   % waveform at end edge of last raster time

% Add time sample 'noise' to test robustness to roundoff error
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;
gNoisy.waveform = g.waveform;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));
gNoisy.first = g.first;
gNoisy.last = g.last;

% Interpolate
rasterOut = 4e-6;
gout = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);

% Plot
tt = (0:(length(gout)-1))*rasterOut;
hold off;
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
hold on;
plot(tt*1e6, gout, 'rx-');
xlabel('time (us)');

% Do it again for a sinusoidal waveform
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;
n = 20;   % number of waveform samples
g.tt = ([1:n]-0.5)*rasterIn;  % times at center of inputer raster intervals
freq = 1000; % Hz 
g.waveform = sin(2*pi*freq*g.tt);
g.first = 0;
g.last = g.waveform(end) + (g.waveform(end)-g.waveform(end-1))/2;

gNoisy.waveform = g.waveform;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));
gNoisy.first = g.first;
gNoisy.last = g.last;

gout = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);
tt = (0:(length(gout)-1))*rasterOut;
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
plot(tt*1e6, gout, 'rx-');

legend('Input waveform', 'Output waveform', 'Input waveform', 'Output waveform'); 
