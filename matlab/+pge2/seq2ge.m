function seq2ge(seqfile, sysGE, pislquant, PNSwt)
% function seq2ge(seqfile, sysGE, pislquant, [PNSwt])
%
% Convert .seq file to a .pge file for the pge2 GE interpreter.
%
% Inputs
%  seqfile    string     Pulseq file
%  sysGE      struct     GE hardware struct, see pge2.opts()
%  pislquant  [1]        number of ADC events at start of scan for receive gain calibration
%  PNSwt      [3]        PNS channel/direction weights. Default: [1 1 1]

if nargin < 4
    PNSwt = [1 1 1];
end

fn = replace(seqfile, {'.seq'}, '');

ceq = seq2ceq([fn '.seq']);

params = pge2.check(ceq, sysGE, 'wt', PNSwt);

pge2.writeceq(ceq, [fn '.pge'], 'pislquant', pislquant, 'params', params);

fprintf('\n');

