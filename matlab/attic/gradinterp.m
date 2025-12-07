function [wavOut, ttOut] = gradinterp(g, rasterIn, rasterOut)
% function [wavOut, ttOut] = gradinterp(g, rasterIn, rasterOut)
%
% Interpolate arbitrary Pulseq gradient to uniform raster.
%
% Inputs:
%  g             Pulseq (arbitrary) gradient struct. Must contain the following fields:
%                g.waveform     gradient waveform samples (a.u.)
%                g.tt           sample times (sec)
%                g.first        waveform value at starting edge of 1st raster time
%                g.last         waveform value at ending edge of last raster time
%  rasterIn      input gradient raster time (sec)
%  rasterOut     output gradient raster time (sec)
%
% Output:
%  gOut          gradient on uniform raster
%  ttOut         output gradient sample times (starts at zero)
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
rasterIn = round(rasterIn/dt)*dt;   % round() is safe. ceil() is not, surprisingly
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

% Output sample times are on regular raster
ttOut = [rasterOut/2:rasterOut:ttIn(end)]';

wavOut = interp1(ttIn, wav, ttOut); %, 'linear', 'extrap');

if any(isnan(wavOut))
    error('gradinterp(): NaN after interpolation');
end

return


function sub_test
% Try a few waveforms to see how the interpolation behaves

rasterIn = 10e-6;   % Siemens gradient raster
rasterOut = 4e-6;   % GE gradient raster

% Add time 'noise' to test robustness to roundoff error
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;

% Define test waveform and sample times.
g.waveform = [1 1 0];   % gradient waveform, a.u.
g.tt = [rasterIn/2 58e-6 20*rasterIn];    % sample times, sec
g.first = g.waveform(1);    % waveform at start edge of first raster time
g.last = g.waveform(end);   % waveform at end edge of last raster time

gNoisy = g;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));

% Interpolate and plot
[gOut ttOut] = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);

hold off;
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
hold on;
plot(ttOut*1e6, gOut, 'rx-');
xlabel('time (us)');
ylabel('Waveform amplitude (a.u.)');

% Do it again for a sinusoidal waveform,
% time samples at center of input raster times
tNoise = 1e-11;   % sec
rasterInNoisy = rasterIn + tNoise;
n = 20;   % number of waveform samples
g.tt = ([1:n]-0.5)*rasterIn;
freq = 2000; % Hz 
g.waveform = sin(2*pi*freq*g.tt);
g.first = 0;
g.last = g.waveform(end) + (g.waveform(end)-g.waveform(end-1))/2;

gNoisy = g;
gNoisy.tt = g.tt + tNoise*randn(1, length(g.tt));

[gOut, ttOut] = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
plot(ttOut*1e6, gOut, 'rx-');

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

[gOut ttOut] = gradinterp(gNoisy, rasterInNoisy, rasterOut + tNoise);
plot(g.tt*1e6, g.waveform, 'bo', 'MarkerSize', 8);
plot(ttOut*1e6, gOut, 'rx-');

legend('Input waveform', 'Output waveform', '', '', '', '');
