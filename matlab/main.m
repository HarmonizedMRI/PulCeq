% create gre2d.seq
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
write2DGRE;

% convert to Ceq struct
ceq = seq2ceq('gre2d.seq');

% write Ceq struct to file. This is what tv7 will load and run.
writeceq(ceq, 'test.tar');

