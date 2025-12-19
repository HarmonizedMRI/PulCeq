function sysGE = opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil, varargin)
% opts - Create struct containing GE MR system hardware parameters
% 
% function sysGE = opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil, varargin)
%
% Inputs:
%   psd_rf_wait     [sec]
%   psd_grd_wait    [sec]
%   b1_max          [Gauss]
%   g_max           [Gauss/cm]
%   slew_max        [Gauss/cm/ms]
%   coil            Gradient coil, see options below. 
%                   Alternatively, specify keyword-arguments 'chronaxie', 'rheobase', 'alpha'
%                   (overrides coil input)
%
% coil       Scanner   Gradient   chronaxie rheobase alpha  gmax  smax
% 'xrmw'     MR750w    XRMW       360d-6    20.0     0.324  33    120
% 'xrm'      MR750     XRM        334d-6    23.4     0.333  50    200
% 'whole'    HDx       TRM WHOLE  370d-6    23.7     0.344  23    77
% 'zoom'     HDx       TRM ZOOM   354d-6    29.1     0.309  40    150
% 'hrmbuhp'  UHP       HRMB       359d-6    26.5     0.370  100   200
% 'hrmw'     Premier   HRMW       642.4d-6  17.9     0.310  70    200
% 'magnus'   MAGNUS    MAGNUS     611d-6    52.2     0.324  300   750
%
% These values are on the scanner in /w/config/Scandbdt.cfg or GRSubsystemHWO.xml
% (e.g., /export/home/mx/host/config/current/GRSubsystemHWO.xml)
% (alpha = EffectivedBdTlength<X,Y,Z>/100)
% See also /w/config/GradientConfig.cfg

% Keyword-argument inputs with defaults. Some were determined empirically using WTools.
% You probably don't want to change these.
arg.GRAD_UPDATE_TIME = 4e-6;    % block duration, and gradient waveform duration, must be an integer multiple of this value
arg.RF_UPDATE_TIME = 2e-6;      % RF pulse duration must be an integer multiple of this value
arg.adc_raster_time = 2e-6;     % ADC dwell/sample time must be an integer multiple of this value
arg.adc_dead_time = 40e-6;      % time required to turn on ADC board
arg.adc_ringdown_time = 0e-6;   % time required to turn off ADC board
arg.rf_dead_time = 72e-6;       % time required to turn on RF amplifier
arg.rf_ringdown_time = 56e-6;   % time required to turn off RF amplifier
arg.segment_dead_time = 12e-6;        % dead time at start of segment
arg.segment_ringdown_time = 105e-6;   % ssi time, plus 4us SSP pulse, plus 1us (probably EOS bit)
arg.gamma = 4.2576e3;                 % Hz/G

arg = vararg_pair(arg, varargin);   % in ../

sysGE.GRAD_UPDATE_TIME = arg.GRAD_UPDATE_TIME;
sysGE.RF_UPDATE_TIME = arg.RF_UPDATE_TIME;
sysGE.adc_raster_time = arg.adc_raster_time;
sysGE.adc_dead_time = arg.adc_dead_time;
sysGE.adc_ringdown_time = arg.adc_ringdown_time;
sysGE.rf_dead_time = arg.rf_dead_time;
sysGE.rf_ringdown_time = arg.rf_ringdown_time;
sysGE.segment_dead_time = arg.segment_dead_time;
sysGE.segment_ringdown_time = arg.segment_ringdown_time;
sysGE.psd_rf_wait = psd_rf_wait;
sysGE.psd_grd_wait = psd_grd_wait;
sysGE.b1_max = b1_max;
sysGE.g_max = g_max;
sysGE.slew_max = slew_max;
sysGE.gamma = arg.gamma;

% PNS model parameters
switch lower(coil)
    case 'xrmw',  chronaxie=360d-6; rheobase=20.0; alpha=0.324;
    case 'xrm',   chronaxie=334d-6; rheobase=23.4; alpha=0.333;
    case 'whole', chronaxie=370d-6; rheobase=23.7; alpha=0.344;
    case 'zoom',  chronaxie=354d-6; rheobase=29.1; alpha=0.309;
    case 'hrmbuhp',  chronaxie=359d-6; rheobase=26.5; alpha=0.370;
    case 'hrmw',  chronaxie=642.4d-6; rheobase=17.9; alpha=0.310;
    case 'magnus', chronaxie=611d-6; rheobase=55.2; alpha=0.324;
    otherwise, error('gradient coil (%s) unkown', coil);
end

sysGE.chronaxie = chronaxie;
sysGE.rheobase = rheobase;
sysGE.alpha = alpha;

assert(abs(sysGE.b1_max) < 1, ...
    'b1_max must be specified in unit of G (typically this limit is 0.25 or a bit smaller)');
