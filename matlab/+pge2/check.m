function params = check(ceq, sysGE, varargin)
%
% Check compatibility of a PulCeq (Ceq) sequence object with the 
% GE scanner specifications in 'sysGE'.
%
% The following are checked:
%  - Sequence block timing
%  - Peak b1 and gradient amplitude/slew
%  - PNS (for one segment at a time)
%
% Inputs
%    ceq
%    sysGE
%
% Options
%    wt     [3]   PNS x/y/z/ channel weights. See pge2.pns().
%    
% To determine if the pge2 interpreter output matches
% the original .seq file, use pge2.validate(...)

arg.wt = [1 1 1];

arg = vararg_pair(arg, varargin);   % in ../

tol = 1e-7;   % timing tolerance. Matches 'eps' in the pge2 EPIC code

% initialize return value
params.b1max = 0;      % max RF amplitude
params.gmax = 0;       % max single-axis gradient amplitude [G/cm]
params.smax = 0;       % max single-axis slew rate in sequence, G/cm/ms
params.hash = DataHash(ceq);

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
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
        v = pge2.checksegment(S, sysGE, 'wt', arg.wt);
    catch ME
        error(sprintf('(segment %d, row %d): %s', i, n, ME.message));
    end

    params.b1max = max(params.b1max, v.b1max);
    params.gmax = max(params.gmax, v.gmax);
    params.smax = max(params.smax, v.smax);

    textprogressbar(n/ceq.nMax*100);

    n = n + ceq.segments(i).nBlocksInSegment;
end
textprogressbar((n-1)/ceq.nMax*100);

textprogressbar(' ok'); 
