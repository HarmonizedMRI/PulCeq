function idb = isdelayblock(block)
% function idb = isdelayblock(block)
% 
% Returns true if block contains no waveforms

idb = block.blockDuration > 0 & isempty(block.rf) & isempty(block.adc) & ...
        isempty(block.gx) & isempty(block.gy) & isempty(block.gz);
