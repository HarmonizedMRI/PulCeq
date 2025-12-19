function writeseq(seq, fov, seqName)
% writeseq - Call seq.checkTiming and seq.write
%
% function writeseq(seq, fov, seqName)

% --- Check sequence timing ---
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% --- Output for execution ---
seq.setDefinition('FOV', fov);
seqName = replace(seqName, {'.seq', '.pge'}, '');
seq.setDefinition('Name', seqName);
seq.write([seqName '.seq']);       % Write to pulseq file


