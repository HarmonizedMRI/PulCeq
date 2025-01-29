function gout = trap4ge(gin, commonRasterTime, sys)
% function gout = trap4ge(gin, commonRasterTime, sys)
%
% Extend trapezoid rise/flat/fall times to be on commonRasterTime boundary.
% This ensures that sample points are on both 10us and 4us
% boundary, so that the interpolation from Siemens (10us) to GE (4us)
% raster time is accurate. This becomes particularly important
% in CAIPI EPI where it is crucial to ensure accurate areas of the y/z blips
% after interpolation.
%
% Inputs
%  gin                struct      Pulseq trapezoid gradient struct, created with mr.makeTrapezoid
%  commonRasterTime   [1]         Should be 20e-6 s (common divisor of 10us and 4us)
%  sys                struct      Pulseq system struct

% Design temporary readout, then round rise/flat/fall times to 20us boundary to ensure
% accurate interpolation from 10us to 4us raster time
% The +mr toolbox doesn't let us specify area here,
% so first set ammplitude to dummy value and then scale as needd.

gout = mr.makeTrapezoid(gin.channel, sys, ...
    'amplitude', gin.amplitude, ...   % dummy value
    'delay', gin.delay, ...
    'riseTime', ceil(gin.riseTime/commonRasterTime)*commonRasterTime, ...
    'flatTime', ceil(gin.flatTime/commonRasterTime)*commonRasterTime, ...
    'fallTime', ceil(gin.fallTime/commonRasterTime)*commonRasterTime);

% scale to preserve area
if abs(gin.area) > 1e-6
    gout.amplitude = gout.amplitude * gin.area/gout.area;
end
    
gout.area = gin.area;

