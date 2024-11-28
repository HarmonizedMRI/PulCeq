system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

if 0
    addpath test
    writetestsequence;   % test.seq

    % convert to Ceq struct
    %ceq = seq2ceq('test.seq');
    ceq = seq2ceq(seq);
else
    ceq = seq2ceq('test.seq');
end

% write Ceq struct to file. This is what pge2 will load and run.
writeceq(ceq, 'test.tar');

