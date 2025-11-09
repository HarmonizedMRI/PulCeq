function checkwaveforms(ceq, sysGE, seq, xmlPath, varargin)
% function checkwaveforms(ceq, sysGE, seq, xmlPath, varargin)
%
% Check agreement between scanner/WTools sequence output 
% and the original Pulseq (.seq) object.
%
% Inputs:
%   ceq       struct       Ceq sequence object, see seq2ceq.m
%   sysGE     struct       System hardware info, see pge2.getsys()

arg.row = [];      % row in .seq file. If not specified, check the whole scan.
arg.plot = false;   % plot every segment instance, continue on keyboard press

arg = vararg_pair(arg, varargin);   % in ../

if ~isempty(arg.row)
    arg.plot = true;
    doNextSegment = false;
else
    arg.row = 1;
    doNextSegment = true;
end

if isempty(arg.row)
    doNextSegment = true;
else
end

axesLinked = false;

% Loop over segments
cnt = 0;   % segment instance counter
n = 1;
while n < ceq.nMax % & cnt < 2
    cnt = cnt + 1;

    % determine block range
    i = ceq.loop(n,1);        % segment index
    n1 = n;
    n2 = n - 1 + ceq.segments(i).nBlocksInSegment;

    if n < arg.row
        n = n2 + 1;
        continue;
    end

    % Pulseq waveforms
    w = seq.waveforms_and_times(true, [n1 n2]);

    % pge2 interpreter waveforms (from WTools)
    d = pge2.read_segment_xml(sprintf('%sscan.xml.%04d', xmlPath, cnt));

    % Ceq object waveforms
    L = ceq.loop(n1:n2, :);
    S = pge2.getsegmentinstance(ceq, i, sysGE, L, 'rotate', true, 'interpolate', true);

    plt.tmin = 0;
    plt.tmax = 0;

    %%
    %% check gradient waveforms 
    %%
    ax = {'gx','gy','gz'};
    for iax = 1:length(ax)
        % pge2 interpreter output (.xml files)
        tt.pge2 = d(iax).time(2:end-3)/1e6 - sysGE.segment_dead_time;     
        %[~,I] = unique(tt.pge2);
        %tt.pge2(I+1) = tt.pge2(I+1) + 1e-12;  % offset any duplicate time points (steps)
        g.pge2 = d(iax).value(2:end-3);

        % Pulseq (.seq) object
        tt.seq = w{iax}(1,:);                   % sec
        g.seq = w{iax}(2,:)/sysGE.gamma/100;    % Gauss/cm

        % Ceq object (after seq2ceq.m conversion)
        tt.ceq = S.(ax{iax}).t - sysGE.segment_dead_time;
        g.ceq = S.(ax{iax}).signal;

        plt.tmin = min(plt.tmin, min(tt.pge2(1), tt.seq(1)));
        plt.tmax = max(plt.tmax, max(tt.pge2(end), tt.seq(end)));

        % max percent difference
        g.seqip = interp1(tt.seq, g.seq, tt.pge2, 'linear', 'extrap');
        e = 100*max(abs(g.seqip-g.pge2))/max(max(abs(g.seqip), 1e-3));

        if e > 3
            fprintf('Gradient waveform mismatch (segment at row %d, %s)\n', n, ax{iax});
            arg.plot = true;
            doNextSegment = false;
        end

        if arg.plot
            subplot(5,1,iax);
            plot(tt.seq, g.seq, 'black-');   % Pulseq (ideal) waveform
            hold on;
            plot(tt.pge2, g.pge2, 'r.-');
            %plot(tt.ceq, g.ceq, 'o');
            hold off
            %legend('pge2', 'ceq', 'seq');
            legend('seq', 'pge2');
            ylabel(sprintf('%s\n(G/cm)', ax{iax}), 'Rotation', 0);
            if iax == 1
                title(sprintf('segment starting at row %d (count = %d)', n, cnt));
            end
        end
    end

    %%
    %% check RF waveform
    %%

    % Pulseq waveform
    tt.seq = w{4}(1,:);                % sec
    rf.seq = w{4}(2,:)/sysGE.gamma;    % Gauss

    % pge2 interpreter output
    tt.rho = d(5).time(2:end-3)/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    [~,I] = unique(tt.rho);
    tt.rho(I+1) = tt.rho(I+1) + 1e-12;  % offset any duplicate time points (steps)
    rho = d(5).value(2:end-3);                     % a.u.
    rho = rho/max(abs(rho)) * max(abs(rf.seq));    % Gauss

    tt.theta = d(6).time(2:end-3)/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    [~,I] = unique(tt.theta);
    tt.theta(I+1) = tt.theta(I+1) + 1e-12;  % offset any duplicate time points (steps)
    theta = d(6).value(2:end-3)/2^23*pi;  % radians
    theta = interp1(tt.theta, theta, tt.rho, 'linear', 'extrap');

    rf.pge2 = rho .* exp(1i*theta);
    tt.pge2 = tt.rho;

    plt.tmin = min(plt.tmin, min(tt.pge2(1), tt.seq(1)));
    plt.tmax = max(plt.tmax, max(tt.pge2(end), tt.seq(end)));

    % max percent difference
    % Note that the pge2 interpreter conjugates the RF 
    rf.seqi = interp1(tt.seq, rf.seq, tt.pge2, 'linear', 'extrap');
    e = 100*max(abs(rf.seqi-conj(rf.pge2)))/max(max(abs(rf.pge2), 1e-6))

    if e > 1
        fprintf('RF waveform mismatch (segment at row %d)\n', n);
        arg.plot = true;
        doNextSegment = false;
        keyboard
    end

    if arg.plot
        subplot(5,1,4); hold off;
        title(sprintf('|RF|, segment %d', cnt));
        plot(tt.pge2, rho, 'r.-'); hold on;
        plot(tt.seq, abs(rf.seq), 'black');
        legend('pge2', 'Pulseq');
        ylabel(sprintf('|RF|\n(Gauss)'), 'Rotation', 0);

        subplot(5,1,5); hold off;
        title(sprintf('∠(conj(RF)), segment %d', cnt));
        plot(tt.theta, -d(6).value(2:end-3)/2^23*pi, 'r.-'); hold on;
        plot(tt.seq, angle(rf.seq), 'black');
        legend('pge2', 'Pulseq');
        ylabel(sprintf('∠RF\n(radians)'), 'Rotation', 0);

        xlabel('time (sec)');
        xlim([plt.tmin plt.tmax]);
        drawnow

        if ~axesLinked
            for sp = 1:5
                subplot(5, 1, sp); 
                ax{sp} = gca;
                grid on;
            end
            linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'x');  % common zoom setting (along time axis) for all tiles
            axesLinked = true;
        end
    end

    if doNextSegment
        if arg.plot
            input('Press Enter key to plot next segment ', "s");
        end
        n = n2 + 1;
    else
        fprintf('Exiting\n');
        return;
    end
end
