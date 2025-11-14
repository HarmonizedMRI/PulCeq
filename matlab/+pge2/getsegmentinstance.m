function S = getsegmentinstance(ceq, i, sysGE, L, varargin)
% function S = getsegmentinstance(ceq, i, sysGE, L)
%
% Inputs:
%   ceq        struct                                    Ceq sequence object
%   i          [1] int                                   Segment ID     
%   sysGE      struct                                    See pge2.opts()
%   L          [nBlocksInSegment size(ceq.loop,2)]       Dynamic scan loop settings (rows from ceq.loop array)
%
% Input options:
%   'plot'           true/false
%   'durationOnly'   true/false     If true, S only contains duration field
%   'rotate'         true/false     If false, display gradients in logical coordinate frame, i.e., 
%                                   before rotating. If true, interpolated gradients are shown.
%   'interpolate'    true/false     Interpolate to 4us, or display corner points (for extended traps)
%
% Output:
%   S               struct containing segment waveforms:
%     S.rf          RF waveform samples (Gauss) and times (sec)
%     S.gx/gy/gz    Gradient waveform samples (Gauss/cm) and times (sec)
%     S.duration    sec. Includes sysGE.dead_time and sysGE.segment_ringdown_time
%     S.pns         PNS stimulation waveform (percent of threshold). Only calculated if 'interpolate' = true.

arg.plot = false;
arg.durationOnly = false;
arg.rotate = false; 
arg.interpolate = false; 

arg = vararg_pair(arg, varargin);   % in ../

if arg.rotate
    arg.interpolate = true;
end

blockIDs = ceq.segments(i).blockIDs;
parentBlocks = ceq.parentBlocks;

assert(length(blockIDs) == size(L,1), 'size(L,1) must be equal to number of blocks in segment');

% Get segment duration and block boundaries
tic = sysGE.segment_dead_time;      % running time counter marking block boundary
for j = 1:length(blockIDs)
    S.tic(j) = tic;
    tic = tic + L(j, 13);
end
S.tic(end+1) = tic;
S.duration = tic + sysGE.segment_ringdown_time;
if arg.durationOnly
    return;
end

% initialize waveforms
S.rf.signal = [];
S.rf.t = [];       
for ax = {'gx','gy','gz'}
    S.(ax{1}).signal = [];
    S.(ax{1}).t = [];
end

% build segment the same way the interpreter does it
tic = sysGE.segment_dead_time;      % running time counter marking block boundary

for j = 1:length(blockIDs)
    p = blockIDs(j);  % parent block index
    assert(p == L(j, 2), 'parent block index in L does not match block specified in segment definition');

    msg1 = sprintf('block %d of %d (parent block %d; block start time %.e s)', j, length(blockIDs), p, tic);

    % variable delay block
    if p == -1  
        % update time counter (block boundary)
        tic = tic + L(j,13);   % sec
        continue;
    end

    % static delay block
    if p == 0
        % update time counter (block boundary)
        tic = tic + L(j,13);   % sec. Should be same as b.blockDuration
        continue;
    end

    b = parentBlocks(p).block;    % parent block

    n = round(b.blockDuration/sysGE.GRAD_UPDATE_TIME);
    if abs(n - b.blockDuration/sysGE.GRAD_UPDATE_TIME) > 1e-7
        throw(MException('block:duration', sprintf('%s: Parent block duration not on sysGE.GRAD_UPDATE_TIME boundary', msg1))); 
    end

    % RF
    if ~isempty(b.rf)
        % get raster time. Assume arbitrary waveform for now. TODO: support ext trap rf
        raster = diff(b.rf.t);
        raster = round(raster(1)/sysGE.RF_UPDATE_TIME)*sysGE.RF_UPDATE_TIME;

        % time samples and waveform
        S.rf.t = [S.rf.t; tic + sysGE.psd_rf_wait + b.rf.delay + b.rf.t]; % + raster/2];
        S.rf.signal = [S.rf.signal; b.rf.signal/max(abs(b.rf.signal))*L(j,3)/sysGE.gamma] * exp(1i*L(j,4));  % complex, Gauss
    end

    % Gradients
    ax = {'gx','gy','gz'};
    grad_amp_indeces = [6 8 10];
    for iax = 1:3
        g = b.(ax{iax});
        if ~isempty(g)
            if strcmp(g.type, 'trap');
                if g.flatTime > eps
                    tt = [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime]';
                    wav = [0; 1; 1; 0]; % normalized amplitude
                else
                    tt = [0 g.riseTime g.riseTime+g.fallTime]';
                    wav = [0; 1; 0]; % normalized amplitude
                end
            else
                % arbitrary gradient or extended trapezoid
                tt = g.tt;
                wav = g.waveform/max(abs(g.waveform));  % normalized amplitude
            end

            wav = L(j,grad_amp_indeces(iax)) / sysGE.gamma / 100 * wav;  % Gauss/cm

            % append to running waveform
            S.(ax{iax}).t = [S.(ax{iax}).t; tic + g.delay + tt];
            S.(ax{iax}).signal = [S.(ax{iax}).signal; wav];
        end
    end

    tic = tic + b.blockDuration;
end

% Remove duplicate gradient samples (extended trapezoids can start/end on block boundary)
[S.gx.t, ia] = unique(S.gx.t);
S.gx.signal = S.gx.signal(ia);
[S.gy.t, ia] = unique(S.gy.t);
S.gy.signal = S.gy.signal(ia);
[S.gz.t, ia] = unique(S.gz.t);
S.gz.signal = S.gz.signal(ia);

if arg.interpolate
    % Interpolate to 4us and rotate gradients,
    % so we can apply the rotation matrix R directly.
    % NB!! The whole segment gets rotated! 

    % Interpolate gradients to sysGE.GRAD_UPDATE_TIME (= 4us)
    [S.gx.t, S.gx.signal] = sub_interp_grad(S.gx.t, S.gx.signal, S.duration, sysGE);
    [S.gy.t, S.gy.signal] = sub_interp_grad(S.gy.t, S.gy.signal, S.duration, sysGE);
    [S.gz.t, S.gz.signal] = sub_interp_grad(S.gz.t, S.gz.signal, S.duration, sysGE);

    if arg.rotate

        % Get rotation matrix R.
        % The rotation matrix is equal to the last non-identity
        % rotation specified in L (if any).
        I = eye(3);
        Iv = I(:);
        for j = length(blockIDs):-1:1
            Rv = L(j, 15:23);   % R in row-major order
            if ~all(round(1e6*Rv) == 1e6*Iv)
                break;
            end
        end

        Rt = reshape(Rv,3,3);  % R transpose
        R = Rt';

        % Apply rotation
        G = R * [S.gx.signal'; S.gy.signal'; S.gz.signal'];
        S.gx.signal = G(1,:)';
        S.gy.signal = G(2,:)';
        S.gz.signal = G(3,:)';
    end
end    

% Calculate PNS waveform
[Spns.gx.t, Spns.gx.signal] = sub_interp_grad(S.gx.t, S.gx.signal, S.duration, sysGE);
[Spns.gy.t, Spns.gy.signal] = sub_interp_grad(S.gy.t, S.gy.signal, S.duration, sysGE);
[Spns.gz.t, Spns.gz.signal] = sub_interp_grad(S.gz.t, S.gz.signal, S.duration, sysGE);
S.pns.t = (1:length(Spns.gx.t))'*sysGE.GRAD_UPDATE_TIME;
Smin = sysGE.rheobase/sysGE.alpha;
G = [Spns.gx.signal'; Spns.gy.signal'; Spns.gz.signal']/100;  % T/m
[S.pns.signal, p] = pge2.pns(Smin, sysGE.chronaxie, G, sysGE.GRAD_UPDATE_TIME, false); 
S.pns.signal = S.pns.signal(:);

if ~arg.plot
    return;
end

% plot
clear ax

subplot(5,1,1);
ax{1} = gca;
plot([0; S.rf.t; S.duration], [0; abs(S.rf.signal); 0], 'black.');
ylabel('RF (Gauss)'); % ylim([0 1.1]);

subplot(5,1,2);
ax{2} = gca;
plot([0; S.rf.t; S.duration], [0; angle(S.rf.signal); 0], 'blue.');
ylabel('RF angle (radians)'); % ylim([0 1.1]);

sp = 3;
cols = 'rgb';

for d = {'gx','gy','gz'}
    subplot(5,1,sp);
    ax{sp} = gca;
    plot([0; S.(d{1}).t; S.duration], [0; S.(d{1}).signal; 0], [cols(sp-2) '.-']);
    ylabel([d ' (G/cm)']);
    ylim(sysGE.g_max*1.05*[-1 1]);
    sp = sp + 1;
end
xlabel('time (sec)');

linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'x');  % common zoom setting (along time axis) for all tiles

return


function [tt, g] = sub_interp_grad(tt, g, dur, sysGE)
    % Interpolate gradients to uniform raster time (4 us)
    % Inputs:
    %  tt      time samples before interpolation, arbitrary points (sec)
    %  g       gradient sampled at tt (a.u.)
    %  dur     total segment duration, including segment dead/ringdown times (sec)
    %  sysGE   pge2 system struct, see opts.m

    dt = sysGE.GRAD_UPDATE_TIME;  % gradient raster time

    t_start = sysGE.segment_dead_time;
    t_end = dur - sysGE.segment_ringdown_time;
    T = t_start + dt/2:dt:t_end;

    if isempty(tt)
        g = zeros(size(T));
    else
        % add 0 samples at start and end of segment so that 'extrap' below works correctly
        tt = [t_start; tt; t_end];
        g = [0; g; 0];

        % remove duplicate samples (e.g., on block boundaries)
        [tt, ia] = unique(tt);
        g = g(ia);

        g = interp1(tt, g, T, 'linear', 'extrap');
    end

    tt = T(:);
    g = g(:);

    return
