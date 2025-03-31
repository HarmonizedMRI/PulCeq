function S = constructvirtualsegment(blockIDs, parentBlocks, sys, plotSegment)
%
% Inputs:
%   blockIDS       ceq.segments(i).blockIDs
%   parentBlocks   ceq.parentBlocks
%   sys            See pge2.getsys()
%  
% Output:
% S               struct containing segment sequencer waveforms:
%   S.rf          RF waveform samples and times (amplitude normalized to 1)
%   S.SSP         A somewhat loose representation of the hardware instruction signals
%                 on the SSP bus. SSP is represented here as a waveform on 4us raster,
%                 with amplitude HI...
%                  - during RF and ADC events, including the dead and ringdown intervals, and
%                  - at the start of pure delay blocks.
%                 The important thing here is that SSP instructions from different
%                 RF/ADC/pure delay events cannot overlap.
%   S.gx/gy/gz    Gradient waveform samples and times (amplitude normalized to 1)

if nargin < 4
    plotSegment = false;
end

% Get segment duration
S.duration = 0;
for j = 1:length(blockIDs)
    p = blockIDs(j);  % parent block index
    if p == 0  % pure delay block
        S.duration = S.duration + 8e-6;
    else
        S.duration = S.duration + parentBlocks(p).block.blockDuration;
        r(j) = rem(parentBlocks(p).block.blockDuration, sys.GRAD_UPDATE_TIME);
    end
end

% Initialize SSP waveform
% SSP set to HI just means that some hardware instruction is being transmitted
% in connection with an RF event, ADC event, or pure delay block.
HI = 1;              
LO = 0;
n = round(S.duration/sys.GRAD_UPDATE_TIME);
S.SSP.signal = LO*ones(n,1);

% Construct segment waveforms
S.rf.signal = [];
S.rf.t = [];       
for ax = {'gx','gy','gz'}
    S.(ax{1}).signal = [];
    S.(ax{1}).t = [];
end

tic = 0;          % running time counter marking block boundary

for j = 1:length(blockIDs)
    p = blockIDs(j);  % parent block index

    msg1 = sprintf('block %d (parent block %d; block start time %.e s)', j, p, tic);

    if p == 0  % pure delay block
        n1 = tic/sys.GRAD_UPDATE_TIME + 1;
        n2 = n1 + 1;
        S.SSP.signal(round(n1:n2)) = [HI; LO];
        tic = tic + 8e-6;
        for ax = {'gx','gy','gz'}
            S.(ax{1}).slew.normalized.peak(j) = 0;
        end
        continue;
    end

    b = parentBlocks(p).block;    % parent block

    n = round(b.blockDuration/sys.GRAD_UPDATE_TIME);
    if abs(n - b.blockDuration/sys.GRAD_UPDATE_TIME) > 1e-7
        throw(MException('block:duration', sprintf('%s: Parent block duration not on sys.GRAD_UPDATE_TIME boundary', msg1))); 
    end

    if ~isempty(b.rf)
        % get raster time. Assume arbitrary waveform for now. TODO: support ext trap rf
        raster = diff(b.rf.t);
        raster = round(raster(1)/sys.RF_UPDATE_TIME)*sys.RF_UPDATE_TIME;
        S.rf.t = [S.rf.t; tic + sys.psd_rf_wait + b.rf.delay + b.rf.t + raster/2];
        S.rf.signal = [S.rf.signal; b.rf.signal/max(abs(b.rf.signal))];

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

    for ax = {'gx','gy','gz'}
        g = b.(ax{1});
        if ~isempty(g)
            if strcmp(g.type, 'trap');
                tt = [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime]';
                wav = [0; 1; 1; 0]; % normalized amplitude
            else
                % arbitrary gradient or extended trapezoid
                tt = g.tt;
                wav = g.waveform/max(abs(g.waveform));    % normalized waveform

                % If j==1, or previous block is a pure delay block, gradient must start near zero
                %max_delta_g_per_sample = sys.slew_max*sys.GRAD_UPDATE_TIME*1e3;
                %if abs(S.(ax{1}).signal(1)) > max_delta_g_per_sample
                %    throw(MException('grad:start', sprintf('%s: Gradients must be (near) zero at start of segment.', msg1)));
                %end

                % If j==length(blockIDs) or next block is a pure delay block, gradient must end near zero
                %if abs(S.(ax{1}).signal(end)) > max_delta_g_per_sample
                %    throw(MException('grad:end', sprintf('%s: Gradients must be (near) zero at end of segment.', msg1)));
                %end
            end

            % append waveform
            S.(ax{1}).t = [S.(ax{1}).t; tic + g.delay + tt];
            S.(ax{1}).signal = [S.(ax{1}).signal; wav];

            % calculate peak normalized slew rate
            slew = diff(wav)./diff(tt);
            S.(ax{1}).slew.normalized.peak(j) = max(abs(slew));
        else
            S.(ax{1}).slew.normalized.peak(j) = 0;
        end
    end

    tic = tic + b.blockDuration;

    % Check for overlapping SSP messages
    if any(S.SSP.signal > 1.5*HI)
        throw(MException('SSP:overlap', sprintf('%s: SSP messages overlap. Try increasing the separation between RF events, ADC events, and pure delay blocks.', msg1)));
    end
end

if ~plotSegment
    return;
end

% plot
subplot(5,1,1);
plot([0; S.rf.t; S.duration], [0; abs(S.rf.signal); 0], 'black-');
ylabel('RF (a.u.)');  ylim([0 1.1]);
subplot(5,1,2);
n = round(S.duration/sys.GRAD_UPDATE_TIME);
plot(((1:n)-0.5)*sys.GRAD_UPDATE_TIME, S.SSP.signal, 'b.');
ylabel('SSP (a.u.)');  ylim([0 1.1]);
sp = 3;
cols = 'rgb';
for ax = {'gx','gy','gz'}
    subplot(5,1,sp);
    plot([0; S.(ax{1}).t; S.duration], [0; S.(ax{1}).signal; 0], [cols(sp-2) '-']);
    ylabel([ax ' (a.u.)']);
    ylim([-1.2 1.2]);
    sp = sp + 1;
end
xlabel('time (sec)');
