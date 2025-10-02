function W = plot(ceq, sys, varargin)
%
% Plot Ceq sequence object
%
% Inputs:
%   ceq     struct       Ceq sequence object, see seq2ceq.m
%   sys     struct       System hardware info, see pge2.getsys()
%
% Input options:
%  timeRange   [1 2]         Requested start and end times (sec). Actual plot will end on a segment boundary.
%  showBLocks  true/false    Draw vertical lines at block boundaries (slow!)
%  logical     true/false    If true, display gradients in logical coordinate frame, i.e., 
%                            before rotating. Corner/sample points are indicated by circles.
%                            If false, gradients are interpolated to 4us and shown as continuous line.

arg.timeRange = [0 ceq.duration];
arg.showBlocks = false; 
arg.logical = false;

arg = vararg_pair(arg, varargin);   % in ../

if arg.timeRange(1) > ceq.duration
    error('Request time exceeds sequence duration');
end

% get waveforms
%W.SSP.signal = [];
W.tic = [];  % block boundary times
tStart = 0;  % start of plot
W.rf.t = 0; W.rf.signal = 0;
W.gx.t = 0; W.gx.signal = 0;
W.gy.t = 0; W.gy.signal = 0;
W.gz.t = 0; W.gz.signal = 0;

n = 1;    % row counter in ceq.loop
tic = 0;  % running timer marking start of segment instance (sec)
while n < ceq.nMax & tic - eps < min(ceq.duration, arg.timeRange(2))
    % get segment index and dynamic (scan loop) information
    i = ceq.loop(n,1);  % segment index
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

    % Get segment instance and add waveforms to running total
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L, 'logical', arg.logical);
    catch ME
        error(sprintf('(n = %d, i = %d): %s\n', n, i, ME.message));
    end

    W.tic = [W.tic(:); tic + S.tic(:)];

    W.rf.t = [W.rf.t; tic + S.rf.t];
    W.gx.t = [W.gx.t; tic + S.gx.t];
    W.gy.t = [W.gy.t; tic + S.gy.t];
    W.gz.t = [W.gz.t; tic + S.gz.t];

    W.rf.signal = [W.rf.signal; S.rf.signal];
    W.gx.signal = [W.gx.signal; S.gx.signal];
    W.gy.signal = [W.gy.signal; S.gy.signal];
    W.gz.signal = [W.gz.signal; S.gz.signal];

    %W.SSP.signal = [W.SSP.signal; S.SSP.signal];

    tic = tic + S.duration;

    n = n + ceq.segments(i).nBlocksInSegment;
end

duration = tic;

% Remove duplicate gradient samples (extended trapezoids can start/end on block boundary)
[W.gx.t, ia] = unique(W.gx.t);
W.gx.signal = W.gx.signal(ia);
[W.gy.t, ia] = unique(W.gy.t);
W.gy.signal = W.gy.signal(ia);
[W.gz.t, ia] = unique(W.gz.t);
W.gz.signal = W.gz.signal(ia);

% plot

subplot(5,1,1);
ax{1} = gca;
%plot([W.rf.t; duration], [abs(W.rf.signal); 0], 'black.');
plot(W.rf.t, abs(W.rf.signal), 'black.');
ylabel('RF (Gauss)'); % ylim([0 1.1]);
if arg.showBlocks
    xtickformat('%.6f');
    sub_plotblockboundary(W.tic, max(abs(W.rf.signal)));
    xticks(W.tic);
    xtickangle(45);
end

subplot(5,1,2);
ax{2} = gca;
plot(W.rf.t, angle(W.rf.signal), 'blue.');
ylabel('RF phase (radians)'); % ylim([0 1.1]);
%n = round(duration/sys.GRAD_UPDATE_TIME);
%plot(tStart + ((1:length(W.SSP.signal))-0.5)*sys.GRAD_UPDATE_TIME, W.SSP.signal, 'b.');
%ylabel('SSP (a.u.)');  ylim([0 1.2*max(W.SSP.signal)]);
if arg.showBlocks
    xtickformat('%.6f');
    %sub_plotblockboundary(W.tic, max(abs(W.SSP.signal)));
    sub_plotblockboundary(W.tic, pi);
    xticks(W.tic);
    xtickangle(45);
end

sp = 3;
cols = 'rgb';
for g = {'gx','gy','gz'}
    subplot(5,1,sp);
    ax{sp} = gca;
    p = plot(W.(g{1}).t, W.(g{1}).signal, [cols(sp-2) '.-']);
    %setDataTipFormat(p, '%.6f');   % does not work :(
    ylabel([g ' (G/cm)']);
    ylim(sys.g_max*1.05*[-1 1]);
    if arg.showBlocks
        xtickformat('%.6f');
        sub_plotblockboundary(W.tic, max(abs(W.(g{1}).signal)));
        xticks(W.tic);
        xtickangle(45);
    end
    sp = sp + 1;
end
xlabel('time (sec)');

linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'x');  % common zoom setting (along time axis) for all tiles

return


function sub_plotblockboundary(T, vs)

    hold on;
    for m = 1:length(T)
        plot([T(m) T(m)], [-vs vs], ':', 'linewidth', 0.2, 'color', 'k');
    end

    return


function setDataTipFormat(h, fmt)
    % setDataTipFormat(h, fmt)
    % h   = graphics object handle(s), e.g. line, scatter, bar
    % fmt = sprintf format string, e.g. '%.6f', '%.4g', etc.

    for k = 1:numel(h)
        obj = h(k);

        % --- Modern method (R2019b+ with DataTipTemplate) ---
        if isprop(obj, 'DataTipTemplate')
            try
                rows = obj.DataTipTemplate.DataTipRows;
                for r = 1:numel(rows)
                    rows(r).ValueFormat = fmt;
                end
            catch
                warning('Could not set DataTipTemplate for object %d', k);
            end

        % --- Fallback for older MATLAB (no DataTipTemplate) ---
        else
            dcm = datacursormode(ancestor(obj,'figure'));
            set(dcm, 'UpdateFcn', @(~,event) localUpdate(event, fmt));
        end
    end

    return

function txt = localUpdate(event, fmt)
    pos = get(event, 'Position');
    txt = cell(1,numel(pos));
    for i = 1:numel(pos)
        txt{i} = sprintf(['Value %d: ' fmt], i, pos(i));
    end
    return

