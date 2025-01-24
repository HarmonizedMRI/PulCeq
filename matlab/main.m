system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% addpath test
% writetestsequence;   % test.seq

% Convert .seq file to a PulCeq (Ceq) object
ceq = seq2ceq('test.seq');

% write Ceq struct to file. This is what pge2 will load and run.
writeceq(ceq, 'test.pge');
