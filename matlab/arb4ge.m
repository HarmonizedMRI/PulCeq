function gout = arb4ge(gin, commonRasterTime, sys)
% function gout = arb4ge(gin, commonRasterTime, sys)
%
% Pad arbitrary (regularly sampled) gradient waveform so the duration,
% and the event delay, is an integer multiple of commonRasterTime.
% This operation will generally NOT preserve total area.
% See also trap4ge.m
%
% Inputs
%  gin                struct      Pulseq arbitrary gradient event, sampled on a regular raster time
%  commonRasterTime   [1]         Should be 20e-6 s for GE (common divisor of 10us and 4us)
%  sys                struct      Pulseq system struct
%
% Output
%  gout       struct     Pulseq arbitrary gradient event

assert(strcmp(gin.type, 'grad'), "Input gradient must be of type 'grad'");

% initialize
gout.type = gin.type;
gout.channel = gin.channel;
gout.delay = ceil(gin.delay/commonRasterTime)*commonRasterTime;
gout.first = gin.first;
gout.last = gin.last;

% get input waveform raster time
rasterIn = diff(gin.tt);
assert(all(rasterIn - rasterIn(1) < 1e-10), 'Input waveform is not sampled on a regular raster time');
rasterIn = round(1e6*rasterIn(1))*1e-6;

% pad waveform as needed
nc = ceil(gin.shape_dur/commonRasterTime);  % number of samples if sampled on commonRasterTime
n1 = gin.shape_dur / rasterIn; 
assert(round(n1) == length(gin.waveform), 'input shape_dur inconsistent with waveform and raster time');
n_pad = round(nc*commonRasterTime/rasterIn - n1);
gout.waveform = [gin.waveform; gin.waveform(end)*ones(n_pad,1)];
gout.tt = -rasterIn/2 + rasterIn*[1:(n1+n_pad)];

gout.shape_dur = rasterIn*length(gout.waveform);
gout.area = rasterIn * sum(gout.waveform);

