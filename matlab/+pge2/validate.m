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
    fprintf('Parent blocks passed timing check\n');
else
    fprintf('Parent blocks failed timing check\n');
end

% Check base (virtual) segments.
ok = true;
for i = 1:ceq.nSegments      % we use 'i' to count segments here and in the EPIC code
    try 
        Sv{i} = pge2.constructvirtualsegment(ceq.segments(i).blockIDs, ceq.parentBlocks, sys);
    catch ME
        ok = false;
        fprintf('Error: Base (virtual) segment %d: %s\n', i, ME.message);
    end
end
if ok
    fprintf('Base segments passed timing check\n');
else
    fprintf('Base segments failed timing check\n');
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
while n < ceq.nMax
    i = ceq.loop(n, 1);   % segment id

    for j = 1:ceq.segments(i).nBlocksInSegment
        % rf amp
        assert(ceq.loop(n,3)/sys.gamma < sys.b1_max + eps, ...
            sprintf('segment %d, block %d: RF amp exceeds limit', i, j));

        % gradient amp
        gamp = ceq.loop(n, 6)/sys.gamma/100;
        assert(gamp < sys.g_max + eps, ...
            sprintf('segment %d, block %d: x gradient amp (%.3f) exceeds limit', i, j, gamp));
        gamp = ceq.loop(n, 8)/sys.gamma/100;
        assert(gamp < sys.g_max + eps, ...
            sprintf('segment %d, block %d: y gradient amp (%.3f) exceeds limit', i, j, gamp));
        gamp = ceq.loop(n, 10)/sys.gamma/100;
        assert(gamp < sys.g_max + eps, ...
            sprintf('segment %d, block %d: z gradient amp (%.3f) exceeds limit', i, j, gamp));

        % gradient slew

        n = n + 1;
    end

    % 
end


