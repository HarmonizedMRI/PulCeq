function sys = getsys(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, gamma, varargin)
%
% Inputs:
%  psd_rf_wait    [sec]
%  psd_grd_wait    [sec]
%  b1_max          [Gauss]
%  g_max          [Gauss/cm]
%  slew_max          [Gauss/cm/ms]
%  gamma            [Hz/Gauss]

% defaults. Some are determined empirically using WTools.
% You probably don't want to change these.
arg.GRAD_UPDATE_TIME = 4e-6;  % block duration, and gradient waveform duration, must be an integer multiple of this value
arg.RF_UPDATE_TIME = 2e-6;    % RF pulse duration must be an integer multiple of this value
arg.adc_raster_time = 2e-6;   % ADC dwell/sample time must be an integer multiple of this value
arg.adc_dead_time = 36e-6;     % time required to turn on ADC board
arg.adc_ringdown_time = 0e-6;  % time required to turn off ADC board
arg.rf_dead_time = 72e-6;       % time required to turn on RF amplifier
arg.rf_ringdown_time = 56e-6;   % time required to turn off RF amplifier

arg = vararg_pair(arg, varargin);   % in ../

sys.GRAD_UPDATE_TIME = arg.GRAD_UPDATE_TIME;
sys.RF_UPDATE_TIME = arg.RF_UPDATE_TIME;
sys.adc_raster_time = arg.adc_raster_time;
sys.adc_dead_time = arg.adc_dead_time;
sys.adc_ringdown_time = arg.adc_ringdown_time;
sys.rf_dead_time = arg.rf_dead_time;
sys.rf_ringdown_time = arg.rf_ringdown_time;
sys.psd_rf_wait = psd_rf_wait;
sys.psd_grd_wait = psd_grd_wait;
sys.b1_max = b1_max;
sys.g_max = g_max;
sys.slew_max = slew_max;
sys.gamma = gamma;

%assert(all([sys.psd_rf_wait sys.psd_grd_wait sys.GRAD_UPDATE_TIME sys.RF_UPDATE_TIME sys.adc_raster_time] < 1e-3), ...
%   'Times must be specified in units of sec');
assert(abs(sys.b1_max) < 1, 'b1_max must be specified in unit of G (typically this limit is 0.25 or a bit smaller)');
