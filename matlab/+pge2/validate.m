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

% Check base (virtual) segments.
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
    fprintf('Base segments PASSED timing check\n');
else
    fprintf('Base segments FAILED timing check\n');
end

% Check scan loop.  TODO
% This is where waveform amplitudes are set, so for each segment instance we must check that:
%  - RF amplitude does not exceed hardware limit
%  - gradient amplitude on each axis does not exceed hardware limit
%  - gradient slew rate on each axis does not exceed hardware limit
%  - gradients are continuous across block boundaries

%a = input('Check scan loop? It might take a while. (y/n) ', "S");
%if ~strcmp(a, 'y') 
%    return;
%end

n = 1;   % block (row) number
textprogressbar('checking scan loop: ');
while n < ceq.nMax 
    i = ceq.loop(n, 1);   % segment id

    for j = 1:ceq.segments(i).nBlocksInSegment
        % check rf amp
        assert(ceq.loop(n,3)/sys.gamma < sys.b1_max + eps, ...
            sprintf('segment %d, block %d: RF amp exceeds limit', i, j));

        % check gradient amplitude and slew
        inds = [6 8 10];  % column indices in loop array containing gradient amplitude
        axes = {'gx', 'gy', 'gz'};
        for d = 1:length(axes)
            ax = axes{d};
            gamp = ceq.loop(n, inds(d))/sys.gamma/100;
            assert(gamp < sys.g_max + eps, ...
                sprintf('segment %d, block %d: %s gradient amp (%.3f) exceeds limit', i, j, ax, gamp));
            slew_max = abs(gamp * S{i}.(ax).slew.normalized.peak(j) * 1e-3);    % G/cm/ms
            assert(slew_max < sys.slew_max + eps, ...
                sprintf('segment %d, block %d: %s gradient slew (%.3f G/cm/ms) exceeds limit', i, j, ax, slew_max));
        end

        % gradient slew  TODO

        % gradient continuity across block boundaries  TODO

        textprogressbar(n/ceq.nMax*100);
        n = n + 1;
    end

    % 
end
textprogressbar(' PASSED'); 


