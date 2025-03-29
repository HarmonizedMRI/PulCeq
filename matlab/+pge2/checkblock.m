function [ok, msg] = checkblock(b, sys)
%
% Check that parent block is compatible with GE hardware
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
    % check peak b1
    b1_max = max(abs(b.rf.signal))/sys.gamma;    % Gauss
    if b1_max > sys.b1_max
        ok = false;
        msg = sprintf('RF peak value (%.3f G) exceeds hardware limit (%.3f G)', b1_max, sys.b1_max);
        return;
    end

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

    
