system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

if 0
addpath test
writetestsequence;   % test.seq

% convert to Ceq struct
%ceq = seq2ceq('test.seq');
ceq = seq2ceq(seq);
end

ceq = seq2ceq('test.seq');

% write Ceq struct to file. This is what tv7 will load and run.
writeceq(ceq, 'test.tar');

