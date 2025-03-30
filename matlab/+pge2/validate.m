function validate(ceq, sys)
%
% Check compatibility of a PulCeq (Ceq) sequence object with the 'Pulseq on GE v2' interpreter.

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    try
        pge2.checkblocktiming(b, sys);
    catch ME
        fprintf('Error: parent block %d: %s\n', p, ME.message);
    end
end

% Check segment timing.
% Here we construct base (virtual) segments and check that:
%  - gradients start/end on zero at the start/end of each segment, respectively
%  - RF events, ADC events, and pure delay blocks are not too close together
for i = 1:ceq.nSegments      % we use 'i' to count segments here and in the EPIC code
    try 
        pge2.constructvirtualsegment(ceq.segments(i).blockIDs, ceq.parentBlocks, sys);
    catch ME
        fprintf('Error: Base (virtual) segment %d: %s\n', i, ME.message);
    end
end

% Check scan loop.  TODO
% This is where waveform amplitudes are set, so for each segment instance we must check that:
%  - RF amplitude does not exceed hardware limit
%  - gradient amplitude on each axis does not exceed hardware limit
%  - gradient slew rate on each axis does not exceed hardware limit
%  - gradients are continuous across block boundaries
a = input('Check scan loop? It might take a while. (y/n) ', "S");
if ~strcmp(a, 'y') 
    return;
end

