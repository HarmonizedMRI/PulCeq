function sys = getsys(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, varargin)

arg.GRAD_UPDATE_TIME = 4e-6;  % block duration, and gradient waveform duration, must be an integer multiple of this value
arg.RF_UPDATE_TIME = 2e-6;    % RF pulse duration must be an integer multiple of this value
arg.adc_raster_time = 2e-6;   % ADC dwell/sample time must be an integer multiple of this value

arg = vararg_pair(arg, varargin);   % in ../

sys.GRAD_UPDATE_TIME = arg.GRAD_UPDATE_TIME;
sys.RF_UPDATE_TIME = arg.RF_UPDATE_TIME;
sys.adc_raster_time = arg.adc_raster_time;
sys.psd_rf_wait = psd_rf_wait;
sys.psd_grd_wait = psd_grd_wait;
sys.b1_max = b1_max;
sys.g_max = g_max;
sys.slew_max = slew_max;

assert(all([sys.psd_rf_wait sys.psd_grd_wait sys.GRAD_UPDATE_TIME sys.RF_UPDATE_TIME sys.adc_raster_time] < 1e-3), ...
   'Times must be specified in units of sec');
assert(abs(sys.b1_max) < 1, 'b1_max must be specified in unit of G (typically this limit is 0.25 or a bit smaller)');
