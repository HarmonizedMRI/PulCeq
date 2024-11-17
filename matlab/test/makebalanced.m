function gbal = makebalanced(g, varargin)
% Add a trapezoid at end of g to make the total area zero.
%
% function gbal = makebalanced(g, varargin)
%
% Inputs:
%   g          1D gradient waveform (G/cm)
% Options:
%   maxSlew    G/cm/ms. Default: 10.
%   maxGrad    G/cm. Default: 5.

% parse inputs
arg.maxSlew = 10;
arg.maxGrad = 5;
arg = vararg_pair(arg, varargin);

maxSlew = arg.maxSlew;
maxGrad = arg.maxGrad;

dt = 4e-6;   % gradient sample/raster time (s)

maxSlew = 0.995*maxSlew;    % so it passes hardware checks
maxGrad = 0.995*maxGrad;

% ramp to zero
dg = -sign(g(end))*maxSlew*dt*1e3;      % G/sample
g = [g(:)' g(end):dg:0];

% add balancing trapezoid
area = sum(g)*dt;    % G/cm*sec
gblip = trapwave2(abs(area), maxGrad, maxSlew, dt*1e3);
gbal = [g(:); -sign(area)*gblip(:)];

return
