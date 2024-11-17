% create gre2d.seq
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

writetestsequence; %

% convert to Ceq struct
%ceq = seq2ceq('gre2d.seq');
ceq = seq2ceq(seq);

% write Ceq struct to file. This is what tv7 will load and run.
writeceq(ceq, 'test.tar');

