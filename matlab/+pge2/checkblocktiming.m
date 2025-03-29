function [ok, msg] = checkblocktiming(b, sys)
% function [ok, msg] = checkblocktiming(b, sys)
%
% Check that parent block timing is compatible with GE hardware.
%
% Inputs
%    b        Pulseq block
%    sys      hardware parameters, see pge2.getsys()

ok = true;
msg = [];

% check block duration
res = b.blockDuration/sys.GRAD_UPDATE_TIME;  % 'resolution' in GE speak is the number of raster/update intervals
if abs(res - round(res)) > 1e-6
    ok = false;
    msg = sprintf('Block duration not on GRAD_UPDATE_TIME (%.0e s) boundary', sys.GRAD_UPDATE_TIME);
    return;
end

% check rf event
if ~isempty(b.rf)
    % check sample times
    dt = diff(b.rf.t);
    res = dt/sys.RF_UPDATE_TIME;
    if any( abs(res - round(res)) > 1e-6 )
        ok = false;
        msg = sprintf('RF sample spacings not on RF_UPDATE_TIME (%.0e s) boundary', sys.RF_UPDATE_TIME);
        return;
    end

    % check delay 
    res = b.rf.delay/sys.GRAD_UPDATE_TIME;
    if abs(res - round(res)) > 1e-6
        ok = false;
        msg = sprintf('RF delay (%.3e) not on GRAD_UPDATE_TIME (%.0e s) boundary', b.rf.delay, sys.GRAD_UPDATE_TIME);
        return;
    end
end

% check gradient events
for ax = {'gx','gy','gz'}
    g = b.(ax{1});
    [ok, msg] = sub_checkgradient(g, sys);
    if ~ok; return; end;
end

return


function [ok, msg] = sub_checkgradient(g, sys)

ok = true;
msg = [];

if isempty(g); return; end;

% check delay 
res = g.delay/sys.GRAD_UPDATE_TIME;
if abs(res - round(res)) > 1e-6
    ok = false;
    msg = sprintf('Gradient delay (%.3e) not on GRAD_UPDATE_TIME (%.0e s) boundary', g.delay, sys.GRAD_UPDATE_TIME);
    return;
end

% check sample times
if strcmp(g.type, 'trap')
    dt = [g.riseTime g.flatTime g.fallTime];  % duration of trapezoid segments
else
    dt = diff(g.t);
end
res = dt/sys.GRAD_UPDATE_TIME;
if any( abs(res - round(res)) > 1e-6 )
    ok = false;
    if strcmp(g.type, 'trap')
        msg = sprintf('Trapezoid rise/flat/fall time not on GRAD_UPDATE_TIME (%.0e s) boundary', sys.GRAD_UPDATE_TIME);
    else
        msg = sprintf('Arbitary gradient/extended trapezoid sample spacing not on GRAD_UPDATE_TIME (%.0e s) boundary', sys.GRAD_UPDATE_TIME);
    end
    return;
end

return

