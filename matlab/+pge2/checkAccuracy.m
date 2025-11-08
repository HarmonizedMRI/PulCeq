function checkAccuracy(ceq, sysGE, seq, xmlPath, plt)
%
% Check accuracy of the Ceq sequence represenation
% and/or the scanner/WTools simulator output (.xml files), 
% against the original Pulseq sequence.

% arg = vararg_pair(arg, varargin);   % in ../

% Loop through segments

n = 1;     % block (row) counter 
cnt = 0;   % segment instance counter

while n < ceq.nMax & cnt < 4
    cnt = cnt + 1;

    % determine block range
    i = ceq.loop(n,1);  % segment index
    n1 = n;
    n2 = n - 1 + ceq.segments(i).nBlocksInSegment;

    % Pulseq waveforms
    w = seq.waveforms_and_times(false, [n1 n2]);
    tt.seq = w{1}(1,:);                    % sec
    gx.seq = w{1}(2,:)/sysGE.gamma/100;   % Gauss/cm

    % Ceq object waveforms
    L = ceq.loop(n1:n2, :);
    S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);
    tt.ceq = S.gx.t - sysGE.segment_dead_time;
    gx.ceq = S.gx.signal;                   

    % pge2 interpreter waveforms (from WTools)
    d = pge2.read_segment_xml(sprintf('%sscan.xml.%04d', xmlPath, cnt));
    tt.pge2 = d(1).time(1:end-3)/1e6 - sysGE.segment_dead_time;     
    gx.pge2 = d(1).value(1:end-3);

    %gx = interp1(d(1).time(1:end-3)/1e6 - sysGE.segment_dead_time, d(1).value(1:end-3), ...
    %    w{1}(1,:), 'linear', 'extrap');

    % plot x gradient
    figure; hold on;
    plot(tt.pge2, gx.pge2, 'rx-');
    %plot(tt.ceq, gx.ceq, 'bo-');
    plot(tt.seq, gx.seq, 'g-');
    title(sprintf('Segment instance %d', cnt));
    %legend('pge2', 'ceq', 'seq');
    legend('pge2', 'seq');

    % Check norm of difference between waveforms (approximate)
    %tt.pge2 = .time([1 2:2:end-1])/1e6;
    %gx.pge2 = d(1).value([1 3:2:end]);
    if 0
    gx.seq = interp1(tt.seq, gx.seq, tt.pge2, 'nearest', 'extrap');
    [cnt norm(gx.seq-gx.pge2)/norm(gx.seq)]
    figure; hold on;
    plot(tt.pge2, gx.seq, 'bo-');
    plot(tt.pge2, gx.pge2, 'rx-');
    end

    n = n2 + 1;
end
