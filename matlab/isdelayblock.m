function idb = isdelayblock(block)
% function idb = isdelayblock(block)
% 
% Returns true if block contains no events, and is not a pure trigger block

idb = block.blockDuration > 0 & ...
      isempty(block.rf) & isempty(block.adc) & ...
      isempty(block.gx) & isempty(block.gy) & isempty(block.gz) & ...
     ~isfield(block, 'trig');

assert(~ (idb == true & isfield(block, 'label')), 'Pure delay blocks cannot contain any labels');
%if idb == true & isfield(block, 'label')
%    keyboard
%end

