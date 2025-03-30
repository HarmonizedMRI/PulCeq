function ok = validate(ceq, sys)
%
% Check compatibility of a 'ceq' sequence object with the 'Pulseq on GE v2' interpreter.

% Check parent block timing.
% Parent blocks are 'virtual' (waveform amplitudes are arbitrary/normalized), so only check
% timing here; waveforms will be checked below for each segment instance in the scan loop.
for p = 1:ceq.nParentBlocks         % we use 'p' to count parent blocks here and in the EPIC code
    b = ceq.parentBlocks(p).block;
    [ok, msg] = pge2.checkblocktiming(b, sys);
    if ~ok
        error(sprintf('%s (parent block %d)', msg, p));
    end
end

% Check segment timing.
% Here we construct base (virtual) segments and check that:
%  - gradients start/end on zero at the start/end of each segment, respectively
%  - RF/ADC events have sufficient gaps before/after to allow for RF/ADC 'turn on'/'turn off' times
for i  = 1:ceq.nSegments      % we use 'i' to count segments here and in the EPIC code
    blockIDs = ceq.segments(i).blockIDs;
    %sub_checkvirtualsegment(ceq, i)
end


% Check scan loop.
% This is where waveform amplitudes are set, so for each segment instance we must check that:
%  - RF amplitude does not exceed hardware limit
%  - gradient amplitude on each axis does not exceed hardware limit
%  - gradient slew rate on each axis does not exceed hardware limit
%  - gradients are continuous across block boundaries

return

