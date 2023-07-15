function issame = compareblocks(seq, b1Events, b2Events, n1, n2)
%
% Compare two Pulseq blocks
%
% Inputs
%   seq       Pulseq object
%   n1, n2    block IDs of blocks to be compared

issame = true;

b1 = seq.getBlock(n1);
b2 = seq.getBlock(n2);

% Compare duration
tol = 1;   % us
if abs(b1.blockDuration - b2.blockDuration) > tol
    issame = false; return;
end

% Are gradients non-unique (same shapes)
for ax = {'gx','gy','gz'}
    if ~comparegradients(b1.(ax{1}), b2.(ax{1}))
        issame = false; return;
    end
end

% Are ADC events same and consistent
if ~compareadc(seq, n1, n2)
    issame = false; return;
end

% Are RF events non-unique (same shapes)
if ~comparerf(b1, b2, b1Events, b2Events, seq)
    issame = false; return;
end

return

function issame = comparerf(b1, b2, b1Events, b2Events, seq)

    if isempty(b1.rf) & isempty(b2.rf)
        issame = true; return;
    end

    if xor(isempty(b1.rf), isempty(b2.rf))
        issame = false; return;
    end

    % Set to equal if mag and phase shape IDs match
    RFid1 = b1Events(2);
    RFid2 = b2Events(2);
    RFevent1 = seq.rfLibrary.data(RFid1);
    RFevent2 = seq.rfLibrary.data(RFid2);
    magShapeID1 = RFevent1.array(3);
    magShapeID2 = RFevent2.array(3);
    phaseShapeID1 = RFevent1.array(4);
    phaseShapeID2 = RFevent2.array(4);

    if magShapeID1 ~= magShapeID2 | ...
       phaseShapeID1 ~= phaseShapeID2
        issame = false;
    else
        issame = true;
    end
return

function issame = comparegradients(g1, g2)

    issame = true;

    if isempty(g1) & isempty(g2)
        issame = true; return;
    end

    if xor(isempty(g1), isempty(g2))
        issame = false; return;
    end

    if ~strcmp(g1.type, g2.type)
        issame = false; return;
    end

    if strcmp(g1.type, 'trap')
        if (g1.riseTime ~= g2.riseTime | ...
            g1.flatTime ~= g2.flatTime | ...
            g1.fallTime ~= g2.fallTime | ...
            g1.delay ~= g2.delay)
            issame = false; return;
        end
    else
        if (g1.shape_id ~= g2.shape_id)
            issame = false; return;
        end
    end

return

function issame = compareadc(seq, n1, n2)

    adc1 = seq.getBlock(n1).adc;
    adc2 = seq.getBlock(n2).adc;

    if isempty(adc1) & isempty(adc2)
        issame = true; return;
    end

    if xor(isempty(adc1), isempty(adc2))
        issame = false; return;
    end

    if (adc1.numSamples == adc2.numSamples | ...
        adc1.dwell == adc2.dwell | ...
        adc1.delay == adc2.delay)
        issame = true;
    else
        error(sprintf('blocks %d, %d: ADC numsamples/dwell/delay must be same for all ADC events', n1, n2));
    end

return

