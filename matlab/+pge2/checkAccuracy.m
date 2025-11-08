function checkAccuracy(ceq, sysGE, seq, xmlPath, plt)
%
% Check accuracy of the Ceq sequence represenation
% and/or the scanner/WTools simulator output (.xml files), 
% against the original Pulseq sequence.

if nargin < 5
    plt = true;
end

% arg = vararg_pair(arg, varargin);   % in ../

% Loop through segments

n = 1;     % block (row) counter 
cnt = 0;   % segment instance counter

while n < ceq.nMax & cnt < 2
    cnt = cnt + 1;

    % determine block range
    i = ceq.loop(n,1);  % segment index
    n1 = n;
    n2 = n - 1 + ceq.segments(i).nBlocksInSegment;

    % Pulseq waveforms
    w = seq.waveforms_and_times(true, [n1 n2]);

    % Ceq object waveforms
    L = ceq.loop(n1:n2, :);
    S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);

    % pge2 interpreter waveforms (from WTools)
    d = pge2.read_segment_xml(sprintf('%sscan.xml.%04d', xmlPath, cnt));

    % check gradients
    ax = {'gx','gy','gz'};
    for iax = 1:length(ax)
        tt.pge2 = d(iax).time(2:end-3)/1e6 - sysGE.segment_dead_time;     
        g.pge2 = d(iax).value(2:end-3);

        tt.seq = w{iax}(1,:);                   % sec
        g.seq = w{iax}(2,:)/sysGE.gamma/100;    % Gauss/cm
        g.seqip = interp1(tt.seq, g.seq, tt.pge2, 'linear', 'extrap');

        tt.ceq = S.(ax{iax}).t - sysGE.segment_dead_time;
        g.ceq = S.(ax{iax}).signal;

        % max percent difference
        e = 100*max(abs(g.seqip-g.pge2))/max(max(abs(g.seqip), 1e-3));

        % Due to interpolation/quantization error, difference isn't exactly zero
        % even if the pge2 output matches the .seq file exactly as intended
        tol = 3;  % so we don't worry about differences below tol
        if e > tol   
            fprintf('difference detected (row %d, segment %d, %s) -- exiting\n', n, i, ax{iax});
            if plt
                % plot gradient
                figure; hold on;
                title(sprintf('%s, segment instance %d', ax{iax}, cnt));
                plot(tt.seq, g.seq, 'black-');   % Pulseq (ideal) waveform
                plot(tt.pge2, g.pge2, 'rx');
                %plot(tt.ceq, g.ceq, 'bo-');
                %legend('pge2', 'ceq', 'seq');
                legend('seq', 'pge2');
            end
    %        return;
        end
    end

    % check RF
    tt.seq = w{4}(1,:);                % sec
    rf.seq = w{4}(2,:)/sysGE.gamma;    % Gauss
    tt.pge2 = d(5).time(2:end-3)/1e6 - sysGE.segment_dead_time;
    rf.rho = d(5).value(2:end-3);      % a.u.
    rf.rho = rf.rho/max(abs(rf.rho)) * max(abs(rf.seq));
    rf.theta = d(6).value(2:end-3);
    figure; hold on
    title(sprintf('rf, segment %d', cnt));
    plot(tt.seq, abs(rf.seq), 'black');
    plot(tt.pge2, rf.rho, 'rx');

    n = n2 + 1;
end
