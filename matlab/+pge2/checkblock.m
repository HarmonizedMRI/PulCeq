function [ok, msg] = checkblock(b, sys)

ok = true;
msg = [];

% check block duration
res = b.blockDuration/sys.GRAD_UPDATE_TIME;  % 'resolution' in GE speak is the number of raster/update intervals
if abs(res - round(res)) > 1e-6
    ok = false;
    msg = sprintf('Block duration not on GRAD_UPDATE_TIME (%.1e s) boundary', sys.GRAD_UPDATE_TIME);
end
