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

% If last sample is not at end edge of last rasterIn interval, add g.last
if mod(g.tt(end), rasterIn)
    ttIn = [ttIn; ttIn(end) + rasterIn - mod(ttIn(end), rasterIn)];
    wav = [wav; g.last];
end

% Output sample times on regular raster.
ttOutLast = ceil(ttIn(end)/rasterOut)*rasterOut;  % extend last sample past end of waveform if necessary
ttOut = [0:rasterOut:ttIn(end)]';

wavOut = interp1(ttIn, wav, ttOut); %, 'linear', 'extrap');

% If last output sample is before ttIn(end),
% add a sample and set value to g.last
if ttOut(end) < ttIn(end)
    wavOut = [wavOut; g.last];
end

if any(isnan(wavOut))
    error('gradinterp(): NaN after interpolation');
end

return


function sub_test

rasterIn = 10e-6;   % Siemens gradient raster
rasterOut = 4e-6;   % GE gradient raster

% Add time 'noise' to test robustness to roundoff error
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;

% Define test waveform and sample times.
g.waveform = [1 1 0];   % gradient waveform, a.u.
g.tt = [rasterIn/2 108e-6 20*rasterIn];    % sample times, sec
g.first = g.waveform(1);    % waveform at start edge of first raster time
g.last = g.waveform(end);   % waveform at end edge of last raster time

gNoisy = g;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));

% Interpolate
gout = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);

% Plot
hold off;
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
hold on;
tt = (0:(length(gout)-1))*rasterOut;
plot(tt*1e6, gout, 'rx-');
xlabel('time (us)');

% Do it again for a sinusoidal waveform,
% time samples at center of input raster times
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;
n = 20;   % number of waveform samples
g.tt = ([1:n]-0.5)*rasterIn;
freq = 1000; % Hz 
g.waveform = sin(2*pi*freq*g.tt);
g.first = 0;
g.last = g.waveform(end) + (g.waveform(end)-g.waveform(end-1))/2;

gNoisy = g;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));

gout = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
tt = (0:(length(gout)-1))*rasterOut;
plot(tt*1e6, gout, 'rx-');

% Do it again for a negative trapezoid
% Duration not on 4us boundary
g.waveform = [-0.2 -0.6 -0.6 0];
g.first = g.waveform(1);
g.last = g.waveform(end);
g.riseTime = 50e-06;
g.flatTime = 30e-06;
g.fallTime = 70e-06;
g.tt = [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime];

gNoisy = g;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));

gout = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
tt = (0:(length(gout)-1))*rasterOut;
plot(tt*1e6, gout, 'rx-');

legend('Input waveform', 'Output waveform', '', '', '', '');
