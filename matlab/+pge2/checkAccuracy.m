function checkAccuracy(ceq, sysGE, xmlPath, seq)
%
% Check accuracy of WTools sequence output (.xml files)
% against the original Pulseq sequence, and against the Ceq object.

% arg = vararg_pair(arg, varargin);   % in ../

% Loop through segments

n = 1;     % row counter in ceq.loop
cnt = 1;   % segment instance counter

while n < ceq.nMax  & cnt < 4

    % determine block range
    i = ceq.loop(n,1);  % segment index
    n1 = n;
    n2 = n - 1 + ceq.segments(i).nBlocksInSegment;

    % Pulseq waveforms
    w = seq.waveforms_and_times(false, [n1 n2]);

    % pge2 waveforms (from WTools)
    d = pge2.read_segment_xml(sprintf('%sscan.xml.%04d', xmlPath, cnt));

    % Ceq object waveforms
    L = ceq.loop(n1:n2, :);
    S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);

    % plot x gradient
    figure; hold on;
    plot(w{1}(1,:), w{1}(2,:)/sysGE.gamma/100, 'b-');  % Pulseq
    plot(d(1).time(1:end-3)/1e6, d(1).value(1:end-3), 'rx');
    plot(S.gx.t, S.gx.signal, 'go');
    title(sprintf('Segment instance %d', cnt));

    n = n2 + 1;
    cnt = cnt + 1;
end
