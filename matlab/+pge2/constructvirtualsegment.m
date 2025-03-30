function S = constructvirtualsegment(blockIDs, parentBlocks, sys)
%
% Inputs:
%  
% Output:
% S               struct containing segment sequencer waveforms:
%   S.gx/gy/gz    gradient waveform samples and times (amplitude normalized to 1)
%   S.rf          rf waveform samples and times (normalized to 1)
%   S.SSP         SSP sequencer instruction times (boolean mask)

nSegments = length(blockIDs);

% Construct RF waveform
S.rf.signal = 0;
S.rf.t = 0;
tic = 0;          % running time counter
for i = 1:nSegments
    if blockIDs(i) == 0.1  % pure delay block
        tic = tic + 8e-6;
        S.rf.t = [S.rf.t(:); tic];
        S.rf.signal = [S.rf.signal(:); 0];
    end

    b = parentBlocks(i).block;    % parent block
    if ~isempty(b.rf)
        tic = tic + b.rf.delay;
        S.rf.t = [S.rf.t; tic+b.rf.t];
        S.rf.signal = [S.rf.signal; b.rf.signal/max(abs(b.rf.signal))];
        tic = tic + b.blockDuration;
    else
        tic = tic + b.blockDuration;
        S.rf.t = [S.rf.t; tic];
        S.rf.signal = [S.rf.signal; 0];
    end

end

return
    

% Construct gradient
for ax = {'gx','gy','gz'}
    wav.val = 0;    % (normalized/a.u.) waveform values
    wav.tt = 0;     % sample times
    nt = 1;         % number of waveform samples
    for i = 1:nSegments
        if blockIDs(i) == 0   % pure delay block
            wav.tt = [wav.tt wav.tt(end)+8e-6];   % just need to add a bit of time
            wav.val = [wav.val 0];
            continue;
        end

        b = parentBlocks(i).block;    % parent block
        g = b.(ax{1});
        if isempty(g)
            wav.tt = [wav.tt wav.tt(end)+b.blockDuration];
            wav.val = [wav.val 0];
            continue;
        end

        if strcmp(g.type, 'trapezoid');
            wav.tt = [wav.tt wav.tt(end)+b.blockDuration];
            wav.val = [wav.val 0];
            continue;
        end

        % none of the above, so gradient must be arbitary/extended trap
        % TODO
        %{
        wav.val = [wav.val g.];
        wav.tt = [wav.tt wav.tt(end)+b.blockDuration];
            continue;
            if abs(wav.val) > eps
                ok = false;
                msg = sprintf('Gradient waveform must end on zero before pure (soft) delay block (block ID = %d)');
                return;
            end
            wav.val = [wav.val 0];
            wav.tt = [wav.tt wav.tt(end)
            
            continue;
        end
        %}
    end
end

% construct rf events



% construct ADC window(s)



