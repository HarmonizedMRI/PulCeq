% create gre2d.seq
addpath ~/github/jfnielsen/TOPPEpsdSourceCode/v6/examples/
write2DGRE;

% convert to Ceq struct
ceq = seq2ceq('gre2d.seq');

% write Ceq struct to file. This is what tv7 will load and run.
writeceq(ceq, 'gre2d.ceq');

