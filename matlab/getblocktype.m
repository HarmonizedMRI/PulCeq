function T = getblocktype(block)
% function T = getblocktype(block)
%
% Returns a vector with flags indicating block type/properties.
%
%  T = [<has rf/gradient/adc> ...
%       <has TRID label> ...
%       <has physio1 trigger> ...
%       <is pure delay block>]

% minDelayBlockDuration = 4e-6 - 100*eps; 

% Are waveforms and/or ADC events defined?
T(1) = isempty(block.rf) & isempty(block.adc) & ...
       isempty(block.gx) & isempty(block.gy) & isempty(block.gz);
T(1) = ~T(1);

% Does the block define a TRID (start of segment) label?
T(2) = false;
if isfield(block, 'label')
    for ii = 1:length(block.label)
        if strcmp(block.label(ii).label, 'TRID')
            T(2) = true;
            break;
        end
    end
end

% Does the block contain a cardiac trigger?
T(3) = false;
if isfield(block, 'trig')
    if strcmp(block.trig.channel, 'physio1')
        T(3) = true;
    end
end

% Is this a pure delay block?
T(4) = ~T(1); % & block.blockDuration >= minDelayBlockDuration;

%assert( ~ (T(1) | block.blockDuration < minDelayBlockDuration), ...
%       sprintf('Duration of block containing rf/gradient/adc events must exceed %es', minDelayBlockDuration) );
