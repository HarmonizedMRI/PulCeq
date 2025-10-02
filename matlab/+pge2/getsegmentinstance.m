function S = getsegmentinstance(ceq, i, sys, L, varargin)
% function S = getsegmentinstance(ceq, i, sys, L)
%
% Inputs:
%   ceq      struct                                    Ceq sequence object
%   i        [1] int                                   Segment ID     
%   sys      struct                                    See pge2.getsys()
%   L        [nBlocksInSegment size(ceq.loop,2)]       Dynamic scan loop settings (rows from ceq.loop array)
%
% Input options:
%   'plot'           true/false
%   'durationOnly'   true/false     If true, S only contains duration field
%   'logical'        true/false     If true, display gradients in logical coordinate frame, i.e., 
%                                   before rotating. Corner/sample points are indicated by circles.
%                                   If false, gradients are interpolated to 4us and shown as continuous line.
%
% Output:
%   S               struct containing segment sequencer waveforms:
%     S.rf          RF waveform samples (Gauss) and times (sec)
%     S.gx/gy/gz    Gradient waveform samples (Gauss/cm) and times (sec)
%     S.SSP         A somewhat loose representation of the hardware instruction signals
%                   on the SSP bus. SSP is represented here as a waveform on 4us raster,
%                   with amplitude HI...
%                    - during RF and ADC events, including the dead and ringdown intervals, and
%                    - at the start of soft delay blocks.
%                   The important thing here is that SSP instructions from different
%                   RF/ADC/soft delay events cannot overlap.
%     
%     S.duration    sec. Includes sys.dead_time and sys.segment_ringdown_time

arg.plot = false;
arg.durationOnly = false;
arg.logical = false; 

arg = vararg_pair(arg, varargin);   % in ../

blockIDs = ceq.segments(i).blockIDs;
parentBlocks = ceq.parentBlocks;

assert(length(blockIDs) == size(L,1), 'size(L,1) must be equal to number of blocks in segment');

% Get segment duration and block boundaries
tic = sys.segment_dead_time;      % running time counter marking block boundary
for j = 1:length(blockIDs)
    S.tic(j) = tic;
    tic = tic + L(j, 13);
end
S.tic(end+1) = tic;
S.duration = tic + sys.segment_ringdown_time;
if arg.durationOnly
    return;
end

% Initialize SSP waveform
% SSP set to HI just means that some hardware instruction is being transmitted
% in connection with an RF event, ADC event, or variable delay block.
HI = 1;              
LO = 0;
n = round(S.duration/sys.GRAD_UPDATE_TIME);
S.SSP.signal = LO*ones(n,1);
S.SSP.signal(1:3) = HI;    % segment dead time

% initialize waveforms
S.rf.signal = [];
S.rf.t = [];       
for ax = {'gx','gy','gz'}
    S.(ax{1}).signal = [];
    S.(ax{1}).t = [];
end

% build segment the same way the interpreter does it
tic = sys.segment_dead_time;      % running time counter marking block boundary

for j = 1:length(blockIDs)
    p = blockIDs(j);  % parent block index
    assert(p == L(j, 2), 'parent block index in L does not match block specified in segment definition');

    msg1 = sprintf('block %d of %d (parent block %d; block start time %.e s)', j, length(blockIDs), p, tic);

    % variable delay block
    if p == -1  
        % variable delay block requires a 4us SSP pulse
        n1 = round(tic/sys.GRAD_UPDATE_TIME) + 1;
        S.SSP.signal(n1) = S.SSP.signal(n1) + HI;

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

    n = round(b.blockDuration/sys.GRAD_UPDATE_TIME);
    if abs(n - b.blockDuration/sys.GRAD_UPDATE_TIME) > 1e-7
        throw(MException('block:duration', sprintf('%s: Parent block duration not on sys.GRAD_UPDATE_TIME boundary', msg1))); 
    end

    % RF
    if ~isempty(b.rf)
        % get raster time. Assume arbitrary waveform for now. TODO: support ext trap rf
        raster = diff(b.rf.t);
        raster = round(raster(1)/sys.RF_UPDATE_TIME)*sys.RF_UPDATE_TIME;

        % time samples and waveform
        S.rf.t = [S.rf.t; tic + sys.psd_rf_wait + b.rf.delay + b.rf.t]; % + raster/2];
        S.rf.signal = [S.rf.signal; b.rf.signal/max(abs(b.rf.signal))*L(j,3)/sys.gamma];  % amplitude in Gauss

        % SSP pulses
        n1 = (tic + sys.psd_rf_wait + b.rf.delay - sys.rf_dead_time)/sys.GRAD_UPDATE_TIME + 1;
        n2 = (tic + sys.psd_rf_wait + b.rf.delay + b.rf.t(end) + raster/2 + sys.rf_ringdown_time)/sys.GRAD_UPDATE_TIME;
        if abs(n1 - abs(n1)) > 1e-7 
            throw(MException('rf:starttime', sprintf('%s: RF waveform must start on a sys.GRAD_UPDATE_TIME boundary', msg1)));
        end
        if abs(n2 - abs(n2)) > 1e-7 
            throw(MException('rf:endtime', sprintf('%s: RF waveform must end on a sys.GRAD_UPDATE_TIME boundary', msg1)));
        end
        if n2 > S.duration/sys.GRAD_UPDATE_TIME 
            throw(MException('rf:endofsegment', sprintf('%s: RF ringdown extends past end of segment', msg1)));
        end
        S.SSP.signal(round(n1:n2)) = S.SSP.signal(round(n1:n2)) + HI;
    end

    % ADC (SSP pulses)
    if ~isempty(b.adc)
        n1 = (tic + sys.psd_grd_wait + b.adc.delay - sys.adc_dead_time)/sys.GRAD_UPDATE_TIME + 1;
        n2 = (tic + sys.psd_grd_wait + b.adc.delay + b.adc.dwell*b.adc.numSamples + sys.adc_ringdown_time)/sys.GRAD_UPDATE_TIME;
        if abs(n1 - abs(n1)) > 1e-7 
            throw(MException('adc:starttime', sprintf('%s: ADC window must start on a sys.GRAD_UPDATE_TIME boundary', msg1)));
        end
        if abs(n2 - abs(n2)) > 1e-7 
            throw(MException('adc:endtime', sprintf('%s: ADC window must end on a sys.GRAD_UPDATE_TIME boundary', msg1)));
        end
        if n2 > S.duration/sys.GRAD_UPDATE_TIME 
            throw(MException('adc:endofsegment', sprintf('%s: ADC ringdown extends past end of segment', msg1)));
        end
        S.SSP.signal(round(n1:n2)) = S.SSP.signal(round(n1:n2)) + HI;
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

            wav = L(j,grad_amp_indeces(iax)) / sys.gamma / 100 * wav;  % Gauss/cm

            % append to running waveform
            S.(ax{iax}).t = [S.(ax{iax}).t; tic + g.delay + tt];
            S.(ax{iax}).signal = [S.(ax{iax}).signal; wav];
        end
    end

    tic = tic + b.blockDuration;

    % Check for overlapping SSP messages
    if any(S.SSP.signal > 1.5*HI)
        warning(sprintf('%s: SSP messages overlap. Try increasing the separation between RF events, ADC events, and soft delay blocks.', msg1));
%        throw(MException('SSP:overlap', sprintf('%s: SSP messages overlap. Try increasing the separation between RF events, ADC events, and soft delay blocks.', msg1)));
    end
end

% Remove duplicate gradient samples (extended trapezoids can start/end on block boundary)
[S.gx.t, ia] = unique(S.gx.t);
S.gx.signal = S.gx.signal(ia);
[S.gy.t, ia] = unique(S.gy.t);
S.gy.signal = S.gy.signal(ia);
[S.gz.t, ia] = unique(S.gz.t);
S.gz.signal = S.gz.signal(ia);

if ~arg.logical
    % Interpolate to 4us and rotate gradients.
    % NB!! The whole segment gets rotated! 

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

    % Interpolate gradients to 4us (sys.GRAD_UPDATE_TIME)
    [S.gx.t, S.gx.signal] = sub_interp_grad(S.gx.t, S.gx.signal, S.duration, sys);
    [S.gy.t, S.gy.signal] = sub_interp_grad(S.gy.t, S.gy.signal, S.duration, sys);
    [S.gz.t, S.gz.signal] = sub_interp_grad(S.gz.t, S.gz.signal, S.duration, sys);
end    

if ~arg.plot
    return;
end

% plot
clear ax

subplot(5,1,1);
ax{1} = gca;
plot([0; S.rf.t; S.duration], [0; abs(S.rf.signal); 0], 'black.');
ylabel('RF (Gauss)'); % ylim([0 1.1]);
%set(gca, 'color', bgColor);  
%set(gca, 'XTick', []);

subplot(5,1,2);
ax{2} = gca;
plot([0; S.rf.t; S.duration], [0; angle(S.rf.signal); 0], 'blue.');
ylabel('RF angle (radians)'); % ylim([0 1.1]);
%n = round(S.duration/sys.GRAD_UPDATE_TIME);
%plot(((1:n)-0.5)*sys.GRAD_UPDATE_TIME, S.SSP.signal, 'b.');
%ylabel('SSP (a.u.)');  ylim([0 1.2]);

sp = 3;
cols = 'rgb';

for d = {'gx','gy','gz'}
    subplot(5,1,sp);
    ax{sp} = gca;
    plot([0; S.(d{1}).t; S.duration], [0; S.(d{1}).signal; 0], [cols(sp-2) '.-']);
    ylabel([d ' (G/cm)']);
    ylim(sys.g_max*1.05*[-1 1]);
    sp = sp + 1;
end
xlabel('time (sec)');

linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'x');  % common zoom setting (along time axis) for all tiles

function [tt, g] = sub_interp_grad(tt, g, dur, sys)
    % Interpolate gradients to 4us raster time
    % Inputs:
    %  tt   time samples, arbitrary points (sec)
    %  g    gradient sampled at tt (a.u.)
    %  dur  segment duration (sec)
    %  sys  pge2 system struct, see getsys.m

    dt = sys.GRAD_UPDATE_TIME;  % gradient raster time
    t_start = sys.segment_dead_time;
    t_end = dur - sys.segment_dead_time - sys.segment_ringdown_time;
    T = t_start + dt/2:dt:t_end;

    if isempty(tt)
        g = zeros(size(T));
    else
        tt = [t_start; tt; t_end];
        g = [0; g; 0];

        [tt, ia] = unique(tt);
        g = g(ia);

        g = interp1(tt, g, T, 'linear', 'extrap');
    end

    tt = T(:);
    g = g(:);

return
