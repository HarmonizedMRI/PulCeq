function S = constructvirtualsegment(blockIDs, parentBlocks, sys)
%
% Inputs:
%  
% Output:
% S               struct containing segment sequencer waveforms:
%   S.gx/gy/gz    gradient waveform samples and times (amplitude normalized to 1)
%   S.rf          rf waveform samples and times (amplitude normalized to 1)
%   S.SSP         hardware instruction signals, represented here as a waveform on GRAD_UPDATE_TIME raster.
%
% Example:
%  >> sys = pge2.getsys(150e-6, 150e-6, 0.25, 5, 20, 4.2576e3);
%  >> S = pge2.constructvirtualsegment(ceq.segments(1).blockIDs, ceq.parentBlocks, sys);

nSegments = length(blockIDs);

% Get segment duration
S.duration = 0;
for i = 1:nSegments
    p = blockIDs(i);  % parent block index
    if p == 0  % pure delay block
        S.duration = S.duration + 8e-6;
    else
        S.duration = S.duration + parentBlocks(p).block.blockDuration;
    end
end

% Initialize SSP waveform
% SSP set to HI just means that some hardware instruction is being transmitted -- here we don't care about the actual value
HI = 1;              
LO = 0;
n = S.duration/sys.GRAD_UPDATE_TIME;
S.SSP.signal = LO*ones(n,1);

% Construct segment waveforms
S.rf.signal = [];
S.rf.t = [];       
for ax = {'gx','gy','gz'}
    S.(ax{1}).signal = [];
    S.(ax{1}).t = [];
end

tic = 0;          % running time counter marking block boundary
for i = 1:nSegments
    p = blockIDs(i);  % parent block index

    if p == 0  % pure delay block
        n1 = tic/sys.GRAD_UPDATE_TIME + 1;
        n2 = n1 + 1;
        S.SSP.signal(round(n1:n2)) = [HI; LO];
        tic = tic + 8e-6;
        continue;
    end

    b = parentBlocks(p).block;    % parent block

    n = round(b.blockDuration/sys.GRAD_UPDATE_TIME);
    if abs(n - b.blockDuration/sys.GRAD_UPDATE_TIME) > 1e-7
        throw(MException('block:duration', sprintf('Parent block %d duration not on sys.GRAD_UPDATE_TIME boundary', p))); 
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
            throw(MException('rf:starttime', 'RF waveform must start on a sys.GRAD_UPDATE_TIME boundary'));
        end
        if abs(n2 - abs(n2)) > 1e-7 
            throw(MException('rf:endtime', 'RF waveform must end on a sys.GRAD_UPDATE_TIME boundary'));
        end
        if n2 > S.duration/sys.GRAD_UPDATE_TIME 
            throw(MException('rf:endofsegment', 'RF ringdown extends past end of segment'));
        end
        S.SSP.signal(round(n1:n2)) = S.SSP.signal(round(n1:n2)) + HI;
    end

    if ~isempty(b.adc)
        n1 = (tic + sys.psd_grd_wait + b.adc.delay - sys.adc_dead_time)/sys.GRAD_UPDATE_TIME + 1;
        n2 = (tic + sys.psd_grd_wait + b.adc.delay + b.adc.dwell*b.adc.numSamples + sys.adc_ringdown_time)/sys.GRAD_UPDATE_TIME;
        if abs(n1 - abs(n1)) > 1e-7 
            throw(MException('adc:starttime', 'ADC window must start on a sys.GRAD_UPDATE_TIME boundary'));
        end
        if abs(n2 - abs(n2)) > 1e-7 
            throw(MException('adc:endtime', 'ADC window must end on a sys.GRAD_UPDATE_TIME boundary'));
        end
        if n2 > S.duration/sys.GRAD_UPDATE_TIME 
            throw(MException('adc:endofsegment', 'ADC ringdown extends past end of segment'));
        end
        S.SSP.signal(round(n1:n2)) = S.SSP.signal(round(n1:n2)) + HI;
    end

    for ax = {'gx','gy','gz'}
        g = b.(ax{1});
        if ~isempty(g)
            if strcmp(g.type, 'trap');
                S.(ax{1}).t = [S.(ax{1}).t; tic + g.delay + [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime]'];
                S.(ax{1}).signal = [S.(ax{1}).signal; 0; 1; 1; 0]; % normalized amplitude
            end
        end
    end

    tic = tic + b.blockDuration;
end

% Check for overlapping SSP messages
if any(S.SSP.signal > 1.5*HI)
    throw(MException('SSP:overlap', 'SSP messages overlap. Try increasing the separation between RF events, ADC events, and pure delay blocks.'));
end

% plot
subplot(5,1,1);
plot([0; S.rf.t; S.duration], [0; abs(S.rf.signal); 0], 'o');
ylabel('RF (a.u.)');
subplot(5,1,2);
n = S.duration/sys.GRAD_UPDATE_TIME;
plot(((1:n)-0.5)*sys.GRAD_UPDATE_TIME, S.SSP.signal, 'o');
ylabel('SSP (a.u.)');
sp = 3;
for ax = {'gx','gy','gz'}
    subplot(5,1,sp);
    sp = sp + 1;
    plot([0; S.(ax{1}).t; S.duration], [0; S.(ax{1}).signal; 0], 'o-');
    ylabel(ax);
    ylim([0 1.2]);
end
