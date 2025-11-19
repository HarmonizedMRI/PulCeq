function pars = check(ceq, sysGE)
%
% Check compatibility of a PulCeq (Ceq) sequence object with the 
% GE scanner specifications in 'sysGE'.
%
% The following are checked:
%  - Sequence block timing
%  - Peak b1 and gradient amplitude/slew
%  - PNS (for one segment at a time)
%
% To determine if the pge2 interpreter output matches
% the original .seq file, use pge2.validate(...)

tol = 1e-7;   % timing tolerance. Matches 'eps' in the pge2 EPIC code

% initialize return value
pars.b1max = 0;      % max RF amplitude
pars.gmax = 0;       % max single-axis gradient amplitude [G/cm]
pars.smax = 0;       % max single-axis slew rate in sequence, G/cm/ms

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
fprintf('Checking parent blocks timing: ');
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    try
        pge2.checkblocktiming(b, sysGE);
    catch ME
        error('Error in parent block %d: %s\n', p, ME.message);
    end
end

% check all segment instances
n = 1;    % row (block) counter in ceq.loop
textprogressbar('pge2.check(): Checking scan loop: ');
while n < ceq.nMax 
    % get segment instance
    i = ceq.loop(n,1);  % segment index
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);  % dynamic info
    try
        S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    % check it
    try
        v = pge2.checksegment(S, sysGE);
    catch ME
        error(sprintf('(segment %d, row %d): %s', i, n, ME.message));
    end

    textprogressbar(n/ceq.nMax*100);

    n = n + ceq.segments(i).nBlocksInSegment;
end
textprogressbar((n-1)/ceq.nMax*100);

textprogressbar(' PASSED'); 
