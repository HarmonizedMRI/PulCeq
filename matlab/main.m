% create gre2d.seq
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% 2D GRE. No extended trapezoids. Arbitrary RF waveform (sinc).
%write2DGRE;

% 3D GRE. No extended trapezoids. Hard pulse specific on two corner points (extended trap).
write3DGRE;

% convert to Ceq struct
%ceq = seq2ceq('gre2d.seq');
ceq = seq2ceq(seq);

% write Ceq struct to file. This is what tv7 will load and run.
writeceq(ceq, 'test.tar');

