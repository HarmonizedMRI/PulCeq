function gout = exttrap4ge(gin, commonRasterTime, sys)
% function gout = exttrap4ge(gin, commonRasterTime, sys)
%
% Pad extended trapezeoid (irrregularly sampled) gradient waveform so the duration
% of each linear section, and the event delay, 
% is an integer multiple of commonRasterTime.
% See also trap4ge.m and arb4ge.m
%
% Inputs
%  gin                struct      Pulseq arbitrary gradient event, sampled on a regular raster time
%  commonRasterTime   [1]         Should be 20e-6 s for GE (common divisor of 10us and 4us)
%  sys                struct      Pulseq system struct
%
% Output
%  gout       struct     Pulseq extended trapezoid gradient event

assert(strcmp(gin.type, 'grad'), "Input gradient must be of type 'grad'");
assert(round(gin.tt(1)*1e6) == 0, 'First corner point must be at time 0');

% initialize
gout.type = gin.type;
gout.channel = gin.channel;
gout.waveform = gin.waveform;
gout.delay = ceil(gin.delay/commonRasterTime)*commonRasterTime;
gout.first = gin.first;
gout.last = gin.last;

% extend each interval to common raster boundary
gout.waveform(1) = gin.waveform(1);
gout.tt(1) = gin.tt(1);
for n = 2:length(gin.waveform)
    dt = gin.tt(n) - gin.tt(n-1);
    gout.tt(n) = gout.tt(n-1) + ceil(dt/commonRasterTime)*commonRasterTime;
end

gout.shape_dur = gout.tt(end);

% calculate area
gout.area = 0;
for n = 2:length(gout.waveform)
    a1 = gout.waveform(n-1);
    a2 = gout.waveform(n);
    dt = gout.tt(n) - gout.tt(n-1);
    gout.area = gout.area + dt*(a2+a1)/2;
end
