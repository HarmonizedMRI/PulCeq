function validate(ceq, sysGE, seq, xmlPath, varargin)
% function validate(ceq, sysGE, seq, xmlPath, varargin)
%
% Check agreement between MR30.2 scanner/WTools sequence output 
% and the original Pulseq (.seq) object.
%
% Inputs:
%   ceq       struct         Ceq sequence object, see seq2ceq.m
%   sysGE     struct         System hardware info, see pge2.getsys()
%   seq       struct         A Pulseq sequence object
%   xmlPath   string or []   Path to folder containing scan.xml.<xxxx> files.
%                            These files are also used by GE's Pulse View sequence plotter.
% 
% Input options:
%   'row'           [1] or 'all'/[]   Check and plot segment starting at this number in .seq file (default: 'all')
%                                     If row is not on a segment boundary, the following segment will be plotted.
%   'plot'          true/FALSE        Plot each segment (continue to next on pressing 'Enter')
%   'threshRFper'   [1]               RF error tolerance (percent rms error). Default: 10.
%   'b1PlotLim'     [1]               RF plot limit (Gauss). Default: sysGE.b1_max 
%
% Usage:
%   1. Call seq2ceq.m to convert .seq file to ceq
%   2. Simulate the .pge file in WTools, or run in MR30.2 VM/scanner, to create xml files
%   3. Call getsys.m to define sysGE for your scanner.
%   4a. Check and plot the first segment:
%       >> pge2.checkwaveforms(ceq, sysGE, seq, xmlPath, 'row', 1, 'plot', true);
%   4b. Check (and plot) each segment one by one until a waveform mismatch is detected:
%      >> pge2.checkwaveforms(ceq, sysGE, seq, xmlPath, 'row', 'all', 'plot', true);
%   4c. Check all segment instances:
%      >> pge2.checkwaveforms(ceq, sysGE, seq, xmlPath);

% Default options
% Re: thresRFper: Some interpolation error is ok; 
% the main failure modes we're after are things like conj/sign change, and gross timing offsets
arg.row = 'all';      
arg.plot = false;   
arg.threshRFper = 50;  
arg.b1PlotLim = sysGE.b1_max;  % Gauss

arg = vararg_pair(arg, varargin);   % in ../

if ischar(arg.row) | isempty(arg.row)
    arg.row = 1;
    doNextSegment = true;
else
    doNextSegment = false;
end

axesLinked = false;

% Loop over segments
teps = 1e-12;
cnt = 0;   % segment instance counter
n = 1;
if ~arg.plot
    textprogressbar('Checking scan loop: ');
end
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
    % and RF/ADC phase offset
    d = pge2.read_segment_xml(sprintf('%sscan.xml.%04d', xmlPath, cnt));
    th = pge2.readthetaregisters(sprintf('%sscan.xml.%04d.ssp', xmlPath, cnt));
    %phaseOffset.pge2 = th(1).theta/2^23*pi;

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
        tt.pge2 = d(iax).time/1e6 - sysGE.segment_dead_time;     
        plt.tmin = min(plt.tmin, min(tt.pge2(1)));
        plt.tmax = max(plt.tmax, max(tt.pge2(end)));
        g.pge2 = d(iax).value;

        % Pulseq (.seq) object
        tt.seq = w{iax}(1,:);                   % sec
        g.seq = w{iax}(2,:)/sysGE.gamma/100;    % Gauss/cm

        % Ceq object (after seq2ceq.m conversion)
        tt.ceq = S.(ax{iax}).t - sysGE.segment_dead_time;
        g.ceq = S.(ax{iax}).signal;
        I = tt.ceq < tt.pge2(1) | tt.ceq > tt.pge2(end);
        tt.ceq(I) = [];
        g.ceq(I) = [];

        % Check difference with seq object.
        % For traps/ext traps, interpreter waveform is piecewise constant 
        % so allow for small differences due to that fact.
        % In addition, increase tolerance to 1.5x (slew*4us) since seq.waveforms_and_times() 
        % may not be entirely accurate (?)

        if length(tt.seq) > 0
            [g.pge2i, I] = sub_robustinterp1(tt.pge2, g.pge2, tt.seq);
            tmp = g.seq(I);  % if I is full/sparse this is either row/column vector :(
            [err, Imaxdiff] = max(abs(g.pge2i(:)-tmp(:)));    % max difference, G/cm
        else
            err = 0;  % no gradient is present on the current axis
        end
        tol = 3 * sysGE.slew_max * sysGE.GRAD_UPDATE_TIME * 1e3;  % max allowed difference per 4us sample 

        if err > tol
            fprintf('%s waveform mismatch (segment at row %d: max diff %.3f G/cm at t = %.3f ms)\n', ax{iax}, n, err, 1e3*tt.pge2(Imaxdiff));
            doNextSegment = false;
        end

        if arg.plot
            subplot(5,1,iax);
            plot(1e3*tt.pge2, g.pge2, 'r.-');
            hold on;
            plot(1e3*tt.seq, g.seq, 'black-');  
            %plot(tt.ceq, g.ceq, 'g.-');
            hold off
            legend('pge2', 'Pulseq');
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
    rf.seq = w{4}(2,:)/sysGE.gamma;    % complex, Gauss

    % pge2 interpreter waveform
    tt.rho = d(5).time/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    plt.tmin = min(plt.tmin, min(tt.rho(1)));
    plt.tmax = max(plt.tmax, max(tt.rho(end)));
    rho = d(5).value;                     % a.u.
    if max(abs(rho)) > 0                  % avoid divide by zero
        rho = rho/max(abs(rho)) * max(abs(rf.seq));    % Gauss
    end

    tt.theta = d(6).time/1e6 - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    theta = d(6).value/2^23*pi;  % + phaseOffset.pge2;  % radians. TODO: add phase and freq offsets
    theta = angle(exp(-1i*theta));  % minus sign since the pge2 interpreter conjugates the phase

    % Ceq object waveform (output of seq2ceq.m)
    tt.ceq = S.rf.t - sysGE.segment_dead_time - sysGE.psd_rf_wait;
    rf.ceq = S.rf.signal;

    % interpolate and compare
    [thetai, I] = sub_robustinterp1(tt.theta, theta, tt.rho);
    rf.pge2 = rho(I) .* exp(1i*thetai);
    tt.pge2 = tt.rho(I);
    dt = sysGE.GRAD_UPDATE_TIME;

    if length(rf.seq) > 0
        [rf.pge2i, I] = sub_robustinterp1(tt.pge2, rf.pge2, tt.seq);
        tmp = rf.seq(I);  % if I is full/sparse this is either row/column vector :(
        if norm(rf.pge2i) > 0
            %err = 100 * rmse(abs(rf.seqi), abs(rf.pge2(I))) / rmse(rf.seqi, 0*rf.seqi);    % percent rmse
            err = 100 * rmse(abs(rf.pge2i), abs(tmp)) / rmse(rf.pge2i, 0*rf.pge2i);    % percent rmse
        else
            err = 0;
        end
    else
        err = 0;
    end

    if err > arg.threshRFper
        fprintf('RF waveform mismatch (%.1f%%; segment at row %d)\n', err, n);
        doNextSegment = false;
    end

    if arg.plot
        subplot(5,1,4); hold off;
        title(sprintf('|RF|, segment %d', cnt));
        plot(1e3*tt.seq, abs(rf.seq), 'black');
        hold on
        plot(1e3*tt.rho, rho, 'r.-'); 
        %plot(tt.ceq, abs(rf.ceq), 'g.-');
        legend('Pulseq', 'pge2'); 
        ylabel(sprintf('|RF|\n(Gauss)'), 'Rotation', 0);

        subplot(5,1,5); hold off;
        title(sprintf('∠RF, segment %d', cnt));
        plot(1e3*tt.seq, angle(rf.seq), 'black');
        hold on
        plot(1e3*tt.theta, theta, 'r.');
        %plot(tt.ceq, angle(rf.ceq), 'g.-');
        legend('Pulseq', 'pge2'); 
        ylabel(sprintf('∠RF\n(radians)'), 'Rotation', 0);

        xlabel('time (ms)');

        % set plot limits
        for sp = 1:5
            subplot(5, 1, sp);
            xlim(1e3*[plt.tmin plt.tmax]);
            switch sp
                case 4
                    ylim([0 arg.b1PlotLim]);
                case 5
                    ylim(1.1* pi * [-1 1]);
                otherwise
                    ylim(1.1*sysGE.g_max*[-1 1]);
            end
        end

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
        textprogressbar(n/ceq.nMax*100);
        n = n2 + 1;
    else
        fprintf('Exiting\n');
        return;
    end
end

textprogressbar((n-1)/ceq.nMax*100);
fprintf('\n');


% Inputs
%   ttout   [n]   
%   ttin    [m]
%   win     [m]
%
% Output
%   wout    size(ttout)
function [wout, Ikeep] = sub_robustinterp1(ttin, win, ttout)

    [nr, nc] = size(ttout);  

    % make times unique
    ttin = ttin(:) + 1e-12*(0:numel(ttin)-1)';
    ttout = ttout(:) + 1e-12*(0:numel(ttout)-1)';

    w = interp1(ttin, win(:), ttout);

    Ikeep = ~isnan(w);
    wout = w(Ikeep);

    if nr > nc
        wout = wout(:);
    else
        wout = wout(:).';
    end

    return
