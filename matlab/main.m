system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

if 0
    addpath test
    writetestsequence;   % test.seq
    ceq = seq2ceq('test.seq');
else
    %ceq = seq2ceq('~/Downloads/se_30sli_2mm_R3_15dir_b1k.seq');
    ceq = seq2ceq('~/Downloads/epidiff_3_shot_ref_2p0mm_30sli.seq');
end


% write Ceq struct to file. This is what pge2 will load and run.
writeceq(ceq, 'test.tar');
