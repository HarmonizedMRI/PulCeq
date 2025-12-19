function ok = checkblocktiming(b, sysGE)
% checkblocktiming - Check if parent block timing is compatible with GE hardware.
%
% function ok = checkblocktiming(b, sysGE)
%
% Inputs
%    b        Pulseq block
%    sysGE      hardware parameters, see pge2.opts()

ok = false;

% check block duration
res = b.blockDuration/sysGE.GRAD_UPDATE_TIME;  % 'resolution' in GE speak is the number of raster/update intervals
if abs(res - round(res)) > 1e-9
    throw(MException('block:duration', sprintf('Block duration not on GRAD_UPDATE_TIME (%.0e s) boundary', sysGE.GRAD_UPDATE_TIME)));
end

% check rf event
if ~isempty(b.rf)
    % check sample times
    dt = diff(b.rf.t);
    res = dt/sysGE.RF_UPDATE_TIME;
    if any( abs(res - round(res)) > 1e-9 )
        throw(MException('rf:times', sprintf('RF sample spacings not on RF_UPDATE_TIME (%.0e s) boundary', sysGE.RF_UPDATE_TIME)));
    end

    % check delay 
    res = b.rf.delay/sysGE.RF_UPDATE_TIME;
    if abs(res - round(res)) > 1e-9
        throw(MException('rf:delay', sprintf('RF delay (%.3e) not on RF_UPDATE_TIME (%.0e s) boundary', b.rf.delay, sysGE.RF_UPDATE_TIME)));
    end
end

% check gradient events
for ax = {'gx','gy','gz'}
    g = b.(ax{1});
    try 
        sub_checkgradient(g, sysGE);
    catch ME
        throw(MException('grad:timing', sprintf('%s, %s', ax{1}, ME.message)));
    end
end

% check adc event
if ~isempty(b.adc)
    res = b.adc.dwell/sysGE.adc_raster_time;
    if abs(res - round(res)) > 1e-9
        throw(MException('adc:dwell', sprintf('ADC dwell time (%.3e) not an integer multiple of adc raster time (%.1e)', b.adc.dwell, sysGE.adc_raster_time)));
    end

    res = b.adc.delay/1e-6;   % ADC delay must be multiple of 1us
    if abs(res - round(res)) > 1e-9
        throw(MException('adc:delay', sprintf('ADC delay (%.3e) not on 1us boundary', b.adc.delay)));
    end
end

ok = true;

return


function sub_checkgradient(g, sysGE)

if isempty(g); return; end;

% check delay 
res = g.delay/sysGE.GRAD_UPDATE_TIME;
if abs(res - round(res)) > 1e-9
    throw(MException('grad:delay', sprintf('Gradient delay (%.3e) not on GRAD_UPDATE_TIME (%.0e s) boundary', g.delay, sysGE.GRAD_UPDATE_TIME)));
end

% check sample times
if strcmp(g.type, 'trap')
    dt = [g.riseTime g.flatTime g.fallTime];  % duration of trapezoid segments
else
    dt = diff(g.tt);
end
res = dt/sysGE.GRAD_UPDATE_TIME;
if any( abs(res - round(res)) > 1e-9 )
    if strcmp(g.type, 'trap')
        msg = sprintf('Trapezoid rise/flat/fall time not on GRAD_UPDATE_TIME (%.0e s) boundary', sysGE.GRAD_UPDATE_TIME);
    else
        msg = sprintf('Gradient sample spacing not on GRAD_UPDATE_TIME (%.0e s) boundary', sysGE.GRAD_UPDATE_TIME);
    end
    throw(MException('grad:times', sprintf('%s', msg)));
end

return

