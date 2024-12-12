system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

if 1
    addpath test
    writetestsequence;   % test.seq
end

ceq = seq2ceq('test.seq');

% write Ceq struct to file. This is what pge2 will load and run.
writeceq(ceq, 'test.tar');

