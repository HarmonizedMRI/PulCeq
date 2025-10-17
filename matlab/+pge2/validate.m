function validate(ceq, sys)
%
% Check compatibility of a PulCeq (Ceq) sequence object with the 
% GE scanner specifications in 'sys'.
%

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
ok = true;
fprintf('Checking parent blocks timing: ');
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    try
        pge2.checkblocktiming(b, sys);
    catch ME
        ok = false;
        fprintf('Error in parent block %d: %s\n', p, ME.message);
    end
end
if ok
    fprintf('PASSED\n');
else
    fprintf('FAILED\n');
    return;
end

% loop through segment instances
n = 1;    % row counter in ceq.loop
textprogressbar('Checking scan loop: ');
while n < ceq.nMax 
    % get segment instance
    i = ceq.loop(n,1);  % segment index
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);  % dynamic info
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L, 'logical', false);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    % Block boundaries must be on sys.GRAD_UPDATE_TIME boundary
    res = S.tic/sys.GRAD_UPDATE_TIME;
    if any(abs(S.tic/sys.GRAD_UPDATE_TIME - round(S.tic/sys.GRAD_UPDATE_TIME)) > 1e-6)
        error(sprintf('segment %d, instance at row %d: block boundaries not on sys.GRAD_UPDATE_TIME boundary', i, n));
    end

    % check peak b1 amplitude
    if ~isempty(S.rf.signal)
        b1max = max(abs(S.rf.signal))/sys.gamma;
        assert(b1max < sys.b1_max + eps, ...
            sprintf('segment %d, instance at row %d: RF amp (%.3 G) exceeds limit (%.3f)', i, n, b1max, sys.b1_max));
    end

    % check peak gradient amplitude and slew rate
    for ax = {'gx','gy','gz'}
        if ~isempty(S.(ax{1}).signal)
            gmax = max(abs(S.(ax{1}).signal));
            assert(gmax < sys.g_max + eps, ...
                sprintf('segment %d, instance at row %d: %s amp (%.2f G/cm) exceeds limit (%.1f)', i, n, ax{1}, gmax, sys.g_max));
            slew = diff(S.(ax{1}).signal)./diff(S.(ax{1}).t)/1000;  % G/cm/ms
            smax = max(abs(slew));
            assert(smax < sys.slew_max + eps, ...
                sprintf('segment %d, instance at row %d: %s slew rate (%.2f G/cm/ms) exceeds limit (%.1f)', i, n, ax{1}, smax, sys.slew_max));
        end
    end

    textprogressbar(n/ceq.nMax*100);

    n = n + ceq.segments(i).nBlocksInSegment;
end
textprogressbar(n/ceq.nMax*100);
textprogressbar(' PASSED'); 


