function W = plot(ceq, sys, varargin)
% function W = plot(ceq, sys, varargin)
%
% Plot Ceq sequence object
%
% Inputs:
%   ceq     struct       Ceq sequence object, see seq2ceq.m
%   sys     struct       System hardware info, see pge2.getsys()
%
% Input options:
%  timeRange       [2]             Requested start and end times (sec). Actual plot will end on a segment boundary.
%  blockRange      [2]             Requested start and end blocks. Actual plot will end on a segment boundary.
%  showBlocks      TRUE/false      Draw vertical lines at block boundaries (default: true)
%  rotate          true/FALSE      If false, display gradients in logical coordinate frame, i.e., 
%                                  before rotating. If true, only interpolated gradients are shown.
%  interpolate     true/FALSE      Interpolate to 4us, or display corner points (for extended traps)
%
% Output: 
%  W               struct containing the plotted waveforms

arg.timeRange = [0 ceq.duration];
arg.blockRange = [1 ceq.nMax];
arg.showBlocks = true; 
arg.rotate = false;
arg.interpolate = false;

arg = vararg_pair(arg, varargin);   % in ../

if arg.rotate
    if ~arg.interpolate
        warning('Gradient waveforms interpolated to 4us');
    end
    arg.interpolate = true;
end

nSubPlots = 6;

if arg.timeRange(1) > ceq.duration
    error('Request time exceeds sequence duration');
end

% get peak waveform amplitudes (for setting display range)
yLim.rf = max(ceq.loop(:,3))/sys.gamma;
yLim.phs = pi;
yLim.gx = max(abs(ceq.loop(:,6)))/sys.gamma/100;   % Gauss/cm
yLim.gy = max(abs(ceq.loop(:,8)))/sys.gamma/100;   % Gauss/cm
yLim.gz = max(abs(ceq.loop(:,10)))/sys.gamma/100;   % Gauss/cm

% Loop through sequence segments and display 
W.tic = [];  % block boundary times
tStart = 0;  % start of plot
W.rf.t = 0; W.rf.signal = 0;
W.gx.t = 0; W.gx.signal = 0;
W.gy.t = 0; W.gy.signal = 0;
W.gz.t = 0; W.gz.signal = 0;
W.pns.t = 0; W.pns.signal = 0;

n = 1;    % row counter in ceq.loop
nFirst = [];
tic = 0;  % running timer marking start of segment instance (sec)

while n < arg.blockRange(2) & tic - eps < min(ceq.duration, arg.timeRange(2))

    % get segment index 
    i = ceq.loop(n,1);  % segment index

    % skip segment if needed
    n1 = n;
    n2 = n - 1 + ceq.segments(i).nBlocksInSegment;
    if n < arg.blockRange(1);
        n = n2 + 1;
        continue;
    end

    if isempty(nFirst)
        nFirst = n1;
    end

    % get dynamic (scan loop) information
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);

    % get segment instance duration
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L, 'durationOnly', true);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    if arg.timeRange(1) > tic + eps
        % update time and skip to next segment instance
        tic = tic + S.duration;
        W.rf.t = tic; W.rf.signal = 0;
        W.gx.t = tic; W.gx.signal = 0;
        W.gy.t = tic; W.gy.signal = 0;
        W.gz.t = tic; W.gz.signal = 0;
        tStart = tic;
        n = n + ceq.segments(i).nBlocksInSegment;
        continue;
    end

    % Get segment instance and add to plot
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L, 'rotate', arg.rotate, 'interpolate', arg.interpolate);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    sub_addSegmentInstanceToPlot(tic, S, arg.showBlocks, sys, yLim, nSubPlots);

    % Add waveforms to running total (for return value -- not used for plotting)
    W.tic = [W.tic(:); tic + S.tic(:)];

    W.rf.t = [W.rf.t; tic + S.rf.t];
    W.gx.t = [W.gx.t; tic + S.gx.t];
    W.gy.t = [W.gy.t; tic + S.gy.t];
    W.gz.t = [W.gz.t; tic + S.gz.t];
    W.pns.t = [W.pns.t; tic + S.pns.t];

    W.rf.signal = [W.rf.signal; S.rf.signal];
    W.gx.signal = [W.gx.signal; S.gx.signal];
    W.gy.signal = [W.gy.signal; S.gy.signal];
    W.gz.signal = [W.gz.signal; S.gz.signal];
    W.pns.signal = [W.pns.signal; S.pns.signal];

    % Update time counter and go to next segment
    tic = tic + S.duration;
    n = n + ceq.segments(i).nBlocksInSegment;
end

% linx axes
for sp = 1:nSubPlots
    subplot(nSubPlots, 1, sp); 
    ax{sp} = gca;
    grid on;
end
linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5} ax{6}], 'x');  % common zoom setting (along time axis) for all tiles

% set misc figure properties
subplot(nSubPlots,1,1);
msg = sprintf('Displaying blocks %d:%d.\n', nFirst, n2);
if ~arg.rotate
    msg = sprintf('%sLogical coordinates -- gradient rotations not shown.\n', msg);
end
if arg.interpolate
    msg = sprintf('%sGradient waveforms interpolated to 4us.\n', msg);
else
    msg = sprintf('%sWaveform samples/corner points are shown as defined in the original Pulseq file.\n', msg);
end
%msg = [msg 'Vertical lines show block/segment boundaries.'];
title(msg);

yticks(ax{2}, [-pi 0 pi]);
yticklabels(ax{2}, {'-π', '0', 'π'});

return


function sub_plotboundary(T, tp)
    % T    [nt]    block boundary locations (time within sequence)
    % tp   string  'block' or 'segment' 

    if nargin < 2
        tp = 'block';
    end

    T = unique(T);

    hold on;
    for m = 1:length(T)
        if strcmp(tp, 'block')
            xline(T(m), ':', 'linewidth', 0.2, 'color', 'k');
        else
            xline(T(m), '--', 'linewidth', 1.2, 'color', 'r');
        end
    end

    return

function sub_addSegmentInstanceToPlot(tic, S, showBlocks, sysGE, yLim, nSubPlots)
    % tic  [1]      time of start of segment
    % S    struct   segment instance

    subplot(nSubPlots,1,1); hold on;
    plot(1e3*(tic + S.rf.t), abs(S.rf.signal), 'black.');
    ylabel('RF (Gauss)');
    ylabel({'|b1|', 'Gauss'}, 'Rotation', 0); 
    ylim(1.1 * yLim.rf * [-1 1]);
    if showBlocks
       sub_plotboundary(1e3*(tic + S.tic));
       sub_plotboundary(1e3*[tic + S.tic(1)], 'segment');
    end
    %x = 1e3*S.tic;
    %y = 0*x;
    %plot(x, y, '.');
    %for i = 1:numel(x)
    %    text(x(i), y(i)-0.05, sprintf('Block %d', i), ...
    %         'HorizontalAlignment','center');
    %end

    %xtickformat('%d');
    %xticks(1:length(S.tic));
    % xtickangle(45);

    subplot(nSubPlots,1,2); hold on;
    plot(1e3*(tic + S.rf.t), angle(S.rf.signal), 'black.');
    ylabel({'\angleb1', 'rad'}, 'Rotation', 0); 
    ylim(1.1 * yLim.phs * [-1 1]);
    if showBlocks
       sub_plotboundary(1e3*(tic + S.tic));
       sub_plotboundary(1e3*[tic + S.tic(1)], 'segment');
    end

    sp = 3;
    cols = 'rgb';
    for g = {'gx','gy','gz'}
        subplot(nSubPlots,1,sp); hold on;
        p = plot(1e3*(tic + S.(g{1}).t), S.(g{1}).signal, [cols(sp-2) '.-']);
        %setDataTipFormat(p, '%.6f');   % does not work :(
        ylabel({g{1}, 'G/cm'}, 'Rotation', 0);
        yl = 1.05 * yLim.(g{1});
        ylim(yl * [-1 1]);
        if showBlocks
            sub_plotboundary(1e3*(tic + S.tic));
            sub_plotboundary(1e3*[tic + S.tic(1)], 'segment');
        end
        sp = sp + 1;
    end

    % PNS waveform
    ax = subplot(nSubPlots,1,6); hold on;
    plot(1e3*(tic + S.pns.t), S.pns.signal, 'r-');
    if showBlocks
        sub_plotboundary(1e3*(tic + S.tic));
        sub_plotboundary(1e3*[tic + S.tic(1)], 'segment');
    end
    ylabel(sprintf('PNS waveform\n%% of threshold'), 'Rotation', 0);
    yticks(ax, [20:20:100]);
%   yticklabels(ax, {'0', '80', 'π'});

    xlabel('time (ms)');

    return
