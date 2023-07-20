function seq2ge(seqarg, sysGE, ofname)

% Convert the .seq file/object to the PulCeq representation
ceq = seq2ceq(seqarg);

% Write to TOPPE files
ceq2ge(ceq, sysGE, ofname);
