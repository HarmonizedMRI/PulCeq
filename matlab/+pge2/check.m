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
pars.pnsmax.row = 1;   % index of first block in segment instance with max PNS
pars.pnsmax.val = 0;   % maximum gradient-combind PNS
pars.pnsmax.xval = 0;
pars.pnsmax.yval = 0;
pars.pnsmax.zval = 0;

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
ok = true;
fprintf('Checking parent blocks timing: ');
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    try
        pge2.checkblocktiming(b, sysGE);
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
ok = true;
n = 1;    % row counter in ceq.loop
textprogressbar('Checking scan loop: ');
while n < ceq.nMax 
    % get segment instance
    i = ceq.loop(n,1);  % segment index
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);  % dynamic info
    try
        S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    % Block boundaries must be on sysGE.GRAD_UPDATE_TIME boundary
    res = S.tic/sysGE.GRAD_UPDATE_TIME;
    if any(abs(S.tic/sysGE.GRAD_UPDATE_TIME - round(S.tic/sysGE.GRAD_UPDATE_TIME)) > tol)
        error(sprintf('segment %d, instance at row %d: block boundaries not on sysGE.GRAD_UPDATE_TIME boundary', i, n));
    end

    % check peak b1 amplitude
    if ~isempty(S.rf.signal)
        b1max = max(abs(S.rf.signal));  % Gauss
        pars.b1max = max(b1max, pars.b1max);
        assert(b1max < sysGE.b1_max + eps, ...
            sprintf('segment %d, instance at row %d: RF amp (%.3 G) exceeds limit (%.3f)', i, n, b1max, sysGE.b1_max));
    end

    % check peak gradient amplitude and slew rate
    for ax = {'gx','gy','gz'}
        if ~isempty(S.(ax{1}).signal)
            gmax = max(abs(S.(ax{1}).signal)); % G/cm
            pars.gmax = max(gmax, pars.gmax);
            assert(gmax < sysGE.g_max + eps, ...
                sprintf('segment %d, instance at row %d: %s amp (%.2f G/cm) exceeds limit (%.1f)', i, n, ax{1}, gmax, sysGE.g_max));
            slew = diff(S.(ax{1}).signal)./diff(S.(ax{1}).t)/1000;  % G/cm/ms
            smax = max(abs(slew));
            pars.smax = max(smax, pars.smax);
            assert(smax < sysGE.slew_max + eps, ...
                sprintf('segment %d, instance at row %d: %s slew rate (%.2f G/cm/ms) exceeds limit (%.1f)', i, n, ax{1}, smax, sysGE.slew_max));
        end
    end

    % check and record peak PNS
    Smin = sysGE.rheobase/sysGE.alpha;
    G = [S.gx.signal'; S.gy.signal'; S.gz.signal']/100;  % T/m
    try
        [pt, p] = pge2.pns(Smin, sysGE.chronaxie, G, sysGE.GRAD_UPDATE_TIME, false); 
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end
    if ok & max(pt) > 100
        I = find(pt == max(pt));
        fprintf('(row %d, t = %.3f ms, segment %d): PNS exceeds first controlled mode (100%%)!!!\n', ...
            n, I(1)*sysGE.GRAD_UPDATE_TIME*1e3, i);
        ok = false;
    end
    if ok & max(pt) > 80
        I = find(pt == max(pt));
        fprintf('(row %d, segment %d, t = %.3f ms): PNS exceeds normal mode (80%%)!\n', ...
            n, i, I(1)*sysGE.GRAD_UPDATE_TIME*1e3);
        ok = false;
    end
    if max(pt) > pars.pnsmax.val
        pars.pnsmax.val = max(pars.pnsmax.val, max(pt));
        pars.pnsmax.row = n;  
        pars.pnsmax.xval = max(pars.pnsmax.xval, max(abs(p(1,:))));
        pars.pnsmax.yval = max(pars.pnsmax.yval, max(abs(p(2,:))));
        pars.pnsmax.zval = max(pars.pnsmax.zval, max(abs(p(3,:))));
    end

    textprogressbar(n/ceq.nMax*100);

    n = n + ceq.segments(i).nBlocksInSegment;
end
textprogressbar((n-1)/ceq.nMax*100);

if ok
    textprogressbar(' PASSED'); 
else
    textprogressbar(' FAILED'); 
end
