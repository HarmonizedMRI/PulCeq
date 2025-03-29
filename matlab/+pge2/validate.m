function simulate(ceq, sys)
%
% Construct sequence (assemble waveforms) and plot.
% This function also checks compatibility with the 'Pulseq on GE v2' interpreter.
%
% Times are in us throughout (inputs, and in code)

% Check that parent block timing is compatible with GE timing requirements
for p = 1:ceq.nParentBlocks
    b = ceq.parentBlocks(p).block;
    [ok, msg] = pge2.checkblocktiming(b, sys);
    if ~ok
        error(sprintf('%s (parent block %d)', msg, p));
    end
end

% Construct base (virtual) segments. 


