function plot(ceq, sys)
%
% Plot Ceq sequence object
%
% Inputs:
%   ceq     struct       Ceq sequence object, see seq2ceq.m
%   sys     struct       System hardware info, see pge2.getsys()

W.rf.t = 0; W.rf.signal = 0;
W.gx.t = 0; W.gx.signal = 0;
W.gy.t = 0; W.gy.signal = 0;
W.gz.t = 0; W.gz.signal = 0;
W.SSP.signal = [];

% get waveforms
n = 1;    % row counter in ceq.loop
tic = 0;  % running timer marking start of segment instances
while n < ceq.nMax
    i = ceq.loop(n,1);  % segment index
    L = ceq.loop(n:(n-1+ceq.segments(i).nBlocksInSegment), :);
    try
        S = pge2.getsegmentinstance(ceq, i, sys, L, false);
    catch ME
        fprintf('Error (n = %d, i = %d): %s\n', n, i, ME.message);
    end

    W.rf.t = [W.rf.t; tic + S.rf.t];
    W.rf.signal = [W.rf.signal; S.rf.signal];
    W.SSP.signal = [W.SSP.signal; S.SSP.signal];
    W.gx.t = [W.gx.t; tic + S.gx.t];
    W.gx.signal = [W.gx.signal; S.gx.signal];
    W.gy.t = [W.gy.t; tic + S.gy.t];
    W.gy.signal = [W.gy.signal; S.gy.signal];
    W.gz.t = [W.gz.t; tic + S.gz.t];
    W.gz.signal = [W.gz.signal; S.gz.signal];

    tic = tic + S.duration;

    W.rf.t = [W.rf.t; tic]; W.rf.signal = [W.rf.signal; 0];
    W.gx.t = [W.gx.t; tic]; W.gx.signal = [W.gx.signal; 0];
    W.gy.t = [W.gy.t; tic]; W.gy.signal = [W.gy.signal; 0];
    W.gz.t = [W.gz.t; tic]; W.gz.signal = [W.gz.signal; 0];

    n = n + ceq.segments(i).nBlocksInSegment;
end

duration = tic;

% plot
subplot(5,1,1);
ax{1} = gca;
plot([0; W.rf.t; duration], [0; abs(W.rf.signal); 0], 'black.');
ylabel('RF (Gauss)'); % ylim([0 1.1]);
%set(gca, 'color', bgColor);  
%set(gca, 'XTick', []);

subplot(5,1,2);
ax{2} = gca;
n = round(duration/sys.GRAD_UPDATE_TIME);
plot(((1:n)-0.5)*sys.GRAD_UPDATE_TIME, W.SSP.signal, 'b.');
ylabel('SSP (a.u.)');  ylim([0 1.2]);

sp = 3;
cols = 'rgb';

for d = {'gx','gy','gz'}
    subplot(5,1,sp);
    ax{sp} = gca;
    try
        plot([0; W.(d{1}).t; duration], [0; W.(d{1}).signal; 0], [cols(sp-2) '.-']);
    catch
        keyboard
    end
    ylabel([d ' (G/cm)']);
    ylim(sys.g_max*1.05*[-1 1]);
    sp = sp + 1;
end
xlabel('time (sec)');

linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'x');  % common zoom setting (along time axis) for all tiles

