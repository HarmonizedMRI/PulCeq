system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

if 0
    addpath test
    writetestsequence;   % test.seq
    ceq = seq2ceq('test.seq');
else
    ceq = seq2ceq('epidiff_3_shot_ref_2p0mm_60sli.seq');
end


% write Ceq struct to file. This is what pge2 will load and run.
writeceq(ceq, 'test.tar');
