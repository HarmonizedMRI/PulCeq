function g = makeTrapezoid(channel, area, gmax, slew, dt, varargin)
% function g = makeTrapezoid(channel, Area, maxGrad, maxSlew, gradRasterTime)
% 
% Make a symmetric Pulseq trapezoid with requested area.
% 
% Inputs
%   Area             Hz/m*s
%   maxGrad          Hz/m
%   maxSlew          Hz/m/s
%   gradRasterTime   s 
%
% Input options
%   delay    s        Default: 0
%
% Output
%   g        struct   Pulseq gradient event of type 'trap'
%
%    ^    <----- flatTime ----->
%    |    +--------------------+  maxGrad
%    |   /                      \
%    |  /                        \
%    | /                          \
%  0 |+---+------------------------+--> time
%     0   τ                      2*τ+flatTime

% Parse inputs
arg.Delay = 0;
arg = vararg_pair(arg, varargin);   % in ../

% Precompute max ramp duration and area
tauMax   = ceil(gmax/slew/dt) * dt;
areaRamp = gmax * tauMax / 2;

% Initialize
g = struct;

% Determine whether we get a triangle or trapezoid
if 2*areaRamp > area                          % triangle
    g.flatTime = 0;
    g.amplitude = sqrt(area * slew / 2);
    g.riseTime = ceil(g.amplitude/slew/dt) * dt;
    g.fallTime = g.riseTime;
else                                          % trapezoid
    g.amplitude = gmax;
    g.riseTime = tauMax;
    g.fallTime = tauMax;
    g.flatTime = ceil((area - 2*areaRamp) / (gmax * dt)) * dt;
end

% Correct amplitude to match exact desired area
areaApprox = g.amplitude * (g.riseTime + g.flatTime);
g.amplitude = g.amplitude * area / areaApprox;

% return struct
g = struct('type', 'trap', ...
    'channel', channel, ...
    'amplitude', g.amplitude, ...
    'riseTime', g.riseTime, ...
    'flatTime', g.flatTime, ...
    'fallTime', g.fallTime, ...
    'area', g.amplitude*(g.riseTime + g.flatTime), ...
    'flatArea', g.amplitude * g.flatTime, ...
    'delay', arg.Delay, ...
    'first', 0, ...
    'last', 0);

