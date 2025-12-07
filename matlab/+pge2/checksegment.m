function val = checksegment(S, sysGE, varargin)
%
% Inputs
%    S            segment instance, see pge2.getsegmentinstance()
%    sysGE        hardware parameters, see pge2.opts()
%
% Options:
%    wt     [3]   PNS x/y/z/ channel weights. See pge2.pns().
%
% Output
%    val.b1max/gmax/smax    Max b1/gradient amp/slew rate in segment
%    val.pns                PNS waveform

arg.wt = [1 1 1];

arg = vararg_pair(arg, varargin);   % in ../

b1max = 0;
gmax = 0;
smax = 0;

tol = 1e-7;   % timing tolerance. Matches 'eps' in the pge2 EPIC code

% Block boundaries must be on sysGE.GRAD_UPDATE_TIME boundary
res = S.tic/sysGE.GRAD_UPDATE_TIME;
if any(abs(S.tic/sysGE.GRAD_UPDATE_TIME - round(S.tic/sysGE.GRAD_UPDATE_TIME)) > tol)
    throw(MException('hardware:timing', sprintf('block boundaries not on sysGE.GRAD_UPDATE_TIME boundary')));
end

% check peak b1 amplitude
if ~isempty(S.rf.signal)
    b1max = max(b1max, max(abs(S.rf.signal)));  % Gauss
    if b1max > sysGE.b1_max
        throw(MException('hardware:peakb1', sprintf('RF amp (%.3 G) exceeds limit (%.3f)', b1max, sysGE.b1_max)));
    end
else
    b1max = 0;
end

% check peak gradient amplitude and slew rate
for ax = {'gx','gy','gz'}
    if ~isempty(S.(ax{1}).signal)
        gmax = max(gmax, max(abs(S.(ax{1}).signal))); % G/cm
        if gmax > sysGE.g_max
            throw(MException('hardware:gmax', sprintf('%s peak amplitude (%.2f G/cm) exceeds limit (%.1f)', ax{1}, gmax, sysGE.g_max)));
            %sprintf('segment %d, instance at row %d: %s amp (%.2f G/cm) exceeds limit (%.1f)', i, n, ax{1}, gmax, sysGE.g_max));
        end
        slew = diff(S.(ax{1}).signal)./diff(S.(ax{1}).t)/1000;  % G/cm/ms
        smax = max(smax, max(abs(slew)));
        if smax > sysGE.slew_max
            throw(MException('hardware:slew', sprintf('%s slew rate (%.2f G/cm/ms) exceeds limit (%.1f)', ax{1}, smax, sysGE.slew_max)));
        end
    end
end

% check PNS
Smin = sysGE.rheobase/sysGE.alpha;
G = [S.gx.signal'; S.gy.signal'; S.gz.signal']/100;  % T/m
try
    [pt, p] = pge2.pns(Smin, sysGE.chronaxie, G, sysGE.GRAD_UPDATE_TIME, 'wt', arg.wt); 
catch ME
    error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
end
if max(pt) > 100
    I = find(pt == max(pt));
    throw(MException('safety:pns', sprintf('(t = %.3f ms) PNS exceeds first controlled mode (100%%%%)!!!\n', ...
            I(1)*sysGE.GRAD_UPDATE_TIME*1e3)));
end
if max(pt) > 80
    I = find(pt == max(pt));
    throw(MException('safety:pns', sprintf('(t = %.3f ms) PNS exceeds normal mode (80%%%%)!\n', ...
            I(1)*sysGE.GRAD_UPDATE_TIME*1e3)));
end

% return values
val = struct('b1max', b1max, 'gmax', gmax, 'smax', smax, 'pns', pt);

