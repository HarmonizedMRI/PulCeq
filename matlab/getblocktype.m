function T = getblocktype(block)
% function T = getblocktype(block)
%
% Returns a vector indicating block type/properties.
%
%  T = [<has rf/gradient/adc> ...              0 or 1
%       <has TRID label> ...                   0 or 1
%       <has physio1 trigger> ...              0 or 1
%       <pure delay block>          ]          0, 1 (static delay), or 2 (variable delay)

% minDelayBlockDuration = 4e-6 - 100*eps; 

% Are waveforms and/or ADC events defined?
T(1) = isempty(block.rf) & isempty(block.adc) & ...
       isempty(block.gx) & isempty(block.gy) & isempty(block.gz);
T = 1*T;   % non-logical array
%T(1) = ~T(1);

% Does the block define a TRID (start of segment) label?
T(2) = 0;
if isfield(block, 'label')
    for ii = 1:length(block.label)
        if strcmp(block.label(ii).label, 'TRID')
            T(2) = 1;
            break;
        end
    end
end

% Does the block contain a cardiac trigger?
T(3) = 0;
if isfield(block, 'trig')
    if strcmp(block.trig.channel, 'physio1')
        T(3) = 1;
    end
end

% Pure delay
if isfield(block, 'softDelay')
    if T(1) ~= 1
        error('Soft delay block cannot contain other events');
    end
    T(4) = 2;
else
    T(4) = T(1);
end

%assert( ~ (T(1) | block.blockDuration < minDelayBlockDuration), ...
%       sprintf('Duration of block containing rf/gradient/adc events must exceed %es', minDelayBlockDuration) );
