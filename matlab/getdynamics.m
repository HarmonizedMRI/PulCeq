function loop = getdynamics(block, segmentID, parentBlockID, parentBlock)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in physical (Pulseq) units.
%
% Also return the gradient energy on each axis in (G/cm)^2*sec 

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
            energy.(ax{1}) = (g.amplitude/GAM/100)^2 / 3 * g.riseTime ...   % (G/cm)^2*sec
                          + (g.amplitude/GAM/100)^2 * g.flatTime ... 
                          + (g.amplitude/GAM/100)^2 / 3 * g.fallTime;
        else
            % need to check polarity (sign) with respect to parent block
            amp.(ax{1}) = max(abs(g.waveform));
            if any(g.waveform.*parentBlock.(ax{1}).waveform < 0)
                amp.(ax{1}) = -amp.(ax{1});
            end
        end
    end
end

if ~isempty(block.adc)
    recphs = block.adc.phaseOffset;
end

loop = [segmentID parentBlockID ...
        rfamp rfphs rffreq ...
        amp.gx amp.gy amp.gz ...
        recphs ...
        block.blockDuration ...
        energy.gx energy.gy energy.gz];
    
