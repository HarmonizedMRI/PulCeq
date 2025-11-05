function sys = getsys(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil, varargin)
%
% Inputs:
%   psd_rf_wait     [sec]
%   psd_grd_wait    [sec]
%   b1_max          [Gauss]
%   g_max           [Gauss/cm]
%   slew_max        [Gauss/cm/ms]
%   coil            Gradient coil. Optionally specify chronaxies, rheobase, alpha
%                   Options are:
%
% Scanner  Gradient coil             chronaxie rheobase alpha  gmax  smax
% MR750w   XRMW            'xrmw'    360d-6    20.0     0.324  33    120
% MR750    XRM             'xrm'     334d-6    23.4     0.333  50    200
% HDx      TRM WHOLE       'whole'   370d-6    23.7     0.344  23    77
% HDx      TRM ZOOM        'zoom'    354d-6    29.1     0.309  40    150
% UHP      HRMB            'hrmb'    359d-6    26.5     0.370  100   200
% Premier  HRMW            'hrmw'    642.4d-6  17.9     0.310  70    200
% Magnus   MAGNUS          'magnus'  611d-6    52.2     0.324  300   750
%


% defaults. Some are determined empirically using WTools.
% You probably don't want to change these.
arg.GRAD_UPDATE_TIME = 4e-6;  % block duration, and gradient waveform duration, must be an integer multiple of this value
arg.RF_UPDATE_TIME = 2e-6;    % RF pulse duration must be an integer multiple of this value
arg.adc_raster_time = 2e-6;   % ADC dwell/sample time must be an integer multiple of this value
arg.adc_dead_time = 40e-6;     % time required to turn on ADC board
arg.adc_ringdown_time = 0e-6;  % time required to turn off ADC board
arg.rf_dead_time = 72e-6;       % time required to turn on RF amplifier
arg.rf_ringdown_time = 56e-6;   % time required to turn off RF amplifier
arg.segment_dead_time = 12e-6; % dead time at start of segment
arg.segment_ringdown_time = 105e-6;   % ssi time, plus 4us SSP pulse, plus 1us (probably EOS bit)
arg.gamma = 4.2576e3;            % Hz/G
arg.coil = [];

arg = vararg_pair(arg, varargin);   % in ../

assert(~isempty(arg.coil), 'A gradient coil must be specified');

sys.GRAD_UPDATE_TIME = arg.GRAD_UPDATE_TIME;
sys.RF_UPDATE_TIME = arg.RF_UPDATE_TIME;
sys.adc_raster_time = arg.adc_raster_time;
sys.adc_dead_time = arg.adc_dead_time;
sys.adc_ringdown_time = arg.adc_ringdown_time;
sys.rf_dead_time = arg.rf_dead_time;
sys.rf_ringdown_time = arg.rf_ringdown_time;
sys.segment_dead_time = arg.segment_dead_time;
sys.segment_ringdown_time = arg.segment_ringdown_time;
sys.psd_rf_wait = psd_rf_wait;
sys.psd_grd_wait = psd_grd_wait;
sys.b1_max = b1_max;
sys.g_max = g_max;
sys.slew_max = slew_max;
sys.gamma = arg.gamma;
sys.coil = arg.coil;
switch lower(coil)
    case 'xrmw',  chronaxie=360d-6; rheobase=20.0; alpha=0.324;
    case 'xrm',   chronaxie=334d-6; rheobase=23.4; alpha=0.333;
    case 'whole', chronaxie=370d-6; rheobase=23.7; alpha=0.344;
    case 'zoom',  chronaxie=354d-6; rheobase=29.1; alpha=0.309;
    case 'hrmb',  chronaxie=359d-6; rheobase=26.5; alpha=0.370;
    case 'hrmw',  chronaxie=642.4d-6; rheobase=17.9; alpha=0.310;
    case 'magnus', chronaxie=611d-6; rheobase=55.2; alpha=0.324;
    otherwise, error('gradient coil (%s) unkown', coil);
end

assert(abs(sys.b1_max) < 1, 'b1_max must be specified in unit of G (typically this limit is 0.25 or a bit smaller)');
