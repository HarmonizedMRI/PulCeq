function checkwaveforms(ceq, sysGE, seq, xmlPath, plt)
% function checkwaveforms(ceq, sysGE, seq, xmlPath, [plt = false])
%
% Check accuracy of the gradient and RF waveforms
% (timing and shapes) in the Ceq sequence representation
% and/or the scanner/WTools simulator output (.xml files), 
% against the original Pulseq sequence.

if nargin < 5
    plt = false;
end

% Loop over segments
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

    %%
    %% check gradient waveforms 
    %%
    ax = {'gx','gy','gz'};
    for iax = 1:length(ax)
        % pge2 interpreter output (.xml files)
        tt.pge2 = d(iax).time(2:end-3)/1e6 - sysGE.segment_dead_time;     
        g.pge2 = d(iax).value(2:end-3);

        % Pulseq (.seq) object
        tt.seq = w{iax}(1,:);                   % sec
        g.seq = w{iax}(2,:)/sysGE.gamma/100;    % Gauss/cm
        g.seqip = interp1(tt.seq, g.seq, tt.pge2, 'linear', 'extrap');

        % Ceq object (after seq2ceq.m conversion)
        tt.ceq = S.(ax{iax}).t - sysGE.segment_dead_time;
        g.ceq = S.(ax{iax}).signal;

        % max percent difference
        e = 100*max(abs(g.seqip-g.pge2))/max(max(abs(g.seqip), 1e-3));

        % Due to interpolation/quantization error, difference isn't exactly zero
        % even if the pge2 output matches the .seq file exactly as intended
        tol = 3;  % so we don't worry about differences below tol
        if e > tol   
            fprintf('Gradient waveform mismatch (row %d, segment %d, %s) -- exiting\n', n, i, ax{iax});
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
            return;
        end
    end

    %%
    %% check RF waveform
    %%

    % Pulseq waveform
    tt.seq = w{4}(1,:);                % sec
    rf.seq = w{4}(2,:)/sysGE.gamma;    % Gauss

    % pge2 interpreter output
    tt.pge2 = d(5).time(2:end-3)/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    rho = d(5).value(2:end-3);                     % a.u.
    rho = rho/max(abs(rho)) * max(abs(rf.seq));    % Gauss

    tt.theta = d(6).time(2:end-3)/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    [tt.theta, I] = unique(tt.theta);
    theta = d(6).value(2:end-3);
    theta = theta/2^23 * pi;         % radians
    theta = theta(I);
    theta = interp1(tt.theta, theta, tt.pge2, 'linear', 'extrap');

    rf.pge2 = rho .* exp(1i*theta);

    % check difference and exit if too big
    % Note that the pge2 interpreter conjugates the RF 
    rf.seqi = interp1(tt.seq, rf.seq, tt.pge2); %, 'linear', 'extrap');
    e = 100*max(abs(rf.seqi-conj(rf.pge2)))/max(max(abs(rf.pge2), 1e-6)); 

    tol = 3;
    if e < tol
        fprintf('RF waveform mismatch (row %d, segment %d) -- exiting\n', n, i);

        figure; hold on

        subplot(121); hold on;
        title(sprintf('|RF|, segment %d', cnt));
        plot(tt.pge2, abs(rf.pge2), 'r.-');
        plot(tt.seq, abs(rf.seq), 'black');
        legend('pge2', 'Pulseq');

        subplot(122); hold on;
        title(sprintf('angle(conj(RF)), segment %d', cnt));
        plot(tt.pge2, angle(conj(rf.pge2)), 'r.-');
        plot(tt.seq, angle(rf.seq), 'black');
        legend('pge2', 'Pulseq');
        return;
    end

    n = n2 + 1;
end
