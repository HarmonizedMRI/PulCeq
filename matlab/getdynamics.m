function loop = getdynamics(block, blockGroupID, parentBlockID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in physical (Pulseq) units.

% defaults
rfamp = 0;
rfphs = 0;
rffreq = 0;
amp.gx = 0;
amp.gy = 0;
amp.gz = 0;
recphs = 0;

if ~isempty(block.rf)
    rfamp = max(abs(block.rf.signal));
    rfphs = block.rf.phaseOffset;
    rffreq = block.rf.freqOffset;
end

for ax = {'gx','gy','gz'}
    g = block.(ax{1});
    if ~isempty(g)
        if strcmp(g.type, 'trap')
            amp.(ax{1}) = block.(ax{1}).amplitude;
        else
            amp.(ax{1}) = max(abs(block.(ax{1}).waveform));
        end
    end
end

if ~isempty(block.adc)
    recphs = block.adc.phaseOffset;
end

loop = [blockGroupID parentBlockID rfamp rfphs rffreq amp.gx amp.gy amp.gz recphs];
    
