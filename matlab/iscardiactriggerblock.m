function itb = iscardiactriggerblock(block)
% function itb = iscardiactriggerblock(block)
% 
% Returns true if block contains only a physio1 (cardiac) trigger event

if ~isfield(block, 'trig')
    itb = false;
    return;
end

itb = strcmp(block.trig.channel, 'physio1') & ...
      isempty(block.rf) & isempty(block.adc) & ...
      isempty(block.gx) & isempty(block.gy)  & isempty(block.gz) ;

assert(~ (itb == true & isfield(block, 'label')), 'A block containing only a cardiac trigger cannot also contain a label');
