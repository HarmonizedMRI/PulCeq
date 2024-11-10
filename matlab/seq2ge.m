function seq2ge(ifn, ofn)
% Convert a Pulseq file (http://pulseq.github.io/) to a binary 
% sequence representation in the Pulserver binary format.
% For more information about Pulserver, see https://infn-mri.github.io/pulserver/
%
% Inputs
%   ifn = Pulseq (.seq) input file name, e.g., 'gre.seq'
%   ofn = binary output file name, e.g., 'gre.bin'
%
% Example:
%   seq2ge('gre2d.seq', 'gre2d.bin')
%   Then put gre2d.bin on your GE scanner and run with the 'pge2' interpreter

ceq = seq2ceq(ifn);
writeceq(ceq, ofn);
