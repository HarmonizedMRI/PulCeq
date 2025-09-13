function validate(ceq, sys)
%
% Check compatibility of a PulCeq (Ceq) sequence object with the 'Pulseq on GE v2' interpreter.

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
ok = true;
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    try
        pge2.checkblocktiming(b, sys);
    catch ME
        ok = false;
        fprintf('Error: parent block %d: %s\n', p, ME.message);
    end
end
if ok
    fprintf('Parent blocks PASSED timing check\n');
else
    fprintf('Parent blocks FAILED timing check\n');
end

% Check virtual (base) segments.
ok = true;
for i = 1:ceq.nSegments      % we use 'i' to count segments here and in the EPIC code
    try 
        S{i} = pge2.constructvirtualsegment(ceq.segments(i).blockIDs, ceq.parentBlocks, sys);
    catch ME
        ok = false;
        fprintf('Error: Base (virtual) segment %d: %s\n', i, ME.message);
    end
end
if ok
    fprintf('Virtual segments PASSED timing check\n');
else
    fprintf('Virtual segments FAILED timing check\n');
end

% Check scan loop
% This is where waveform amplitudes are set, so for each segment instance we must check that:
%  - RF amplitude does not exceed hardware limit
%  - gradient amplitude and slew rate on each axis do not exceed hardware limits
%  - gradients are continuous across block boundaries
%  - gradients start and end on zero at beginning and end of each segment

%a = input('Check scan loop? It might take a while. (y/n) ', "S");
%if ~strcmp(a, 'y') 
%    return;
%end

n = 1;   % block (row) number
textprogressbar('Checking scan loop: ');
while n < ceq.nMax 
    % get segment index and dynamic (scan loop) information
    i = ceq.loop(n,1);  % segment index
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);

    % Get segment instance
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L);
    catch ME
        fprintf('Error (n = %d, i = %d): %s\n', n, i, ME.message);
    end

    % Gradients must start and end on zero
    for ax = {'gx','gy','gz'}
        assert(abs(S.(ax{1}).signal(1)) < 2*eps, ...
            sprintf('segment %d: %s must be zero at start of segment'));
        assert(abs(S.(ax{1}).signal(end)) < 2*eps, ...
            sprintf('segment %d: %s must be zero at end of segment'));
    end

    % peak b1 amplitude
    assert(max(abs(S.rf.signal))/sys.gamma < sys.b1_max + eps, ...
        sprintf('segment %d: RF amp exceeds limit', i));

    % peak gradient amplitude and slew
    for ax = {'gx','gy','gz'}
        gmax = max(abs(S.(ax{1}).signal));
        assert(gmax < sys.g_max + eps, ...
            sprintf('segment %d: %s amp (%.2f G/cm) exceeds limit (%.1f)', i, ax{1}, gmax, sys.g_max));
        slew = diff(S.(ax{1}).signal)./diff(S.(ax{1}).t)/1000;
        smax = max(abs(slew));
        assert(smax < sys.slew_max + eps, ...
            sprintf('segment %d: %s slew rate (%.2f G/cm/ms) exceeds limit (%.1f)', i, ax{1}, smax, sys.slew_max));
    end

    textprogressbar(n/ceq.nMax*100);

    n = n + ceq.segments(i).nBlocksInSegment;
end
textprogressbar(' PASSED'); 


