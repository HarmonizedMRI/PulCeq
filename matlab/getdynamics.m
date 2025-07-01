function loop = getdynamics(block, segmentID, parentBlockID, physioTrigger, parentBlock)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in physical (Pulseq) units.
%
% Also return the gradient energy on each axis in (G/cm)^2*sec 
%
% This function just inserts identity rotation matrix --
% the calling function is responsible for updating R as needed.

% defaults
rfamp = 0;
rfphs = 0;
rffreq = 0;
amp.gx = 0;
amp.gy = 0;
amp.gz = 0;
recphs = 0;
energy.gx = 0;
energy.gy = 0;
energy.gz = 0;

R = eye(3);

GAM = 4257.6;   % Hz/Gauss

if ~isempty(block.rf)
    rfamp = max(abs(block.rf.signal));
    rfphs = block.rf.phaseOffset;
    rffreq = block.rf.freqOffset;
end

for ax = {'gx','gy','gz'}
    g = block.(ax{1});
    if ~isempty(g)
        if strcmp(g.type, 'trap')
            amp.(ax{1}) = g.amplitude;
            energy.(ax{1}) = (g.amplitude)^2 / 3 * g.riseTime ...   % (Hz/m)^2*sec
                          + (g.amplitude)^2 * g.flatTime ... 
                          + (g.amplitude)^2 / 3 * g.fallTime;
        else  % arbitrary gradient or extended trapezoid
            energy.(ax{1}) = sum((g.waveform(1:end-1)).^2 .* diff(g.tt));
            % energy.(ax{1}) = sum((g.waveform(1:end-1)/GAM/100).^2 .* diff(g.tt));

            amp.(ax{1}) = max(abs(g.waveform)); % correct magnitude but next we check the sign

            % Need to check polarity (sign) with respect to parent block.
            % In the interpreter, the waveform shape that is loaded
            % is pb.g.waveform/max(abs(pb.g.waveform))
            w1 = parentBlock.(ax{1}).waveform;  % Pulseq gradient event, in physical units
            w2 = w1 / max(abs(w1)); % this is the (normalized) shape that is loaded into waveform memory in interpreter
            %mx = max(g.waveform);
            %mn = min(g.waveform);
            %pbmx = max(parentBlock.(ax{1}).waveform);
            %pbmn = min(parentBlock.(ax{1}).waveform);
            %amp.(ax{1}) = max(abs(g.waveform));
            %if any(g.waveform.*parentBlock.(ax{1}).waveform < 0)
            if w2 .* g.waveform < 1
                amp.(ax{1}) = -amp.(ax{1});
            end
        end
    end
end

if ~isempty(block.adc)
    recphs = block.adc.phaseOffset;
    % save ADC frequency for ADC events
    rffreq = block.adc.freqOffset;
end

loop = [segmentID parentBlockID ...
        rfamp rfphs rffreq ...
        amp.gx energy.gx ...
        amp.gy energy.gy ...
        amp.gz energy.gz ...
        recphs ...
        block.blockDuration ...
        physioTrigger ...
        R(:)'];
    
