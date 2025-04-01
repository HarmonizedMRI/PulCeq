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

% Check scan loop
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
textprogressbar('Checking scan loop: ');
while n < ceq.nMax 
    i = ceq.loop(n, 1);   % segment id

    nbis = ceq.segments(i).nBlocksInSegment;

    assert(all(ceq.loop(n:(n+nbis-1),3)/sys.gamma < sys.b1_max + eps), ...
        sprintf('segment %d: RF amp exceeds limit', i));

    % check gradient amplitude and slew
    inds = [6 8 10];  % column indices in loop array containing gradient amplitude
    axes = {'gx', 'gy', 'gz'};
    for d = 1:length(axes)
        ax = axes{d};
        gamp = ceq.loop(n:(n+nbis-1), inds(d))/sys.gamma/100;
        J = find(gamp > sys.g_max);
        if ~isempty(J)
            for ll = 1:length(J)
                j = J(ll);
                fprintf('(block %d) segment %d, block %d: %s gradient amp (%.3f G/cm) exceeds limit\n', n+j-1, i, j, ax, gamp(j));
            end
        end
        slew_max = abs(gamp' .* S{i}.(ax).slew.normalized.peak * 1e-3);    % G/cm/ms
        J = find(slew_max > sys.slew_max);
        if ~isempty(J)
            for ll = 1:length(J)
                j = J(ll);
                fprintf('(block %d) segment %d, block %d: %s gradient slew (%.3f G/cm/ms) exceeds limit\n', n+j-1, i, j, ax, slew_max(j));
            end
        end
    end

    % check gradient continuity across block boundaries  TODO
    n = n + nbis;

    textprogressbar(n/ceq.nMax*100);
end
textprogressbar(' PASSED'); 


