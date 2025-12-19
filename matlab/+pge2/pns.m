function [pt, p] = pns(Smin, c, g, dt, varargin)
% pns - Calculate PNS using nerve impulse response model in IEC 60601-2-33:2022
%
% function [pt, p] = pns(Smin, c, g, dt, ...)
%
% Quick demo:
%   >> pns('test');
%
% Inputs
%   Smin   [1]     stimulation threshold (rheobase) for constant slew of infinite duration
%                  Note: In GE lingo, Smin = r/alpha, r = rheobase, alpha = effective coil length
%   c      [1]     chronaxie (sec). Nerve impulse response time constant
%   g      [3 n]   x/y/z gradient waveform (T/m), for uniform (raster) sampling
%   dt     [1]     gradient raster (sample) time (sec)
%
% Input options
%   wt     [3]     x/y/z/ channel weights. Default: [1 1 1].
%                  From IEC 60601-2-33:2022 section (12): wt = [0.8 1.0 0.7]
%   plt    true/false   plot 
%
% Output
%   pt      [1 n]  channel-combined PNS waveform (% of stimulation threshold)
%   p       [3 n]  PNS waveform on each gradient channel

arg.wt = [1 1 1];
arg.plt = false;

arg = vararg_pair(arg, varargin);   % in ../

if ischar(Smin)
    if strcmp(Smin, 'test')
        pt = sub_test();
        return;
    end
end

assert(size(g,1) < size(g,2), 'g must be size 3xn');
assert(size(g,1) == 3, 'g must be size 3xn');

n = size(g,2);

% Contribution of a slew impulse at time 0 to PNS at time tau
% IEC 60601-2-33:2022 Eq. AA.21
tau = dt:dt:(20*c);  % no need to make much longer than chronaxie
f = dt/Smin * c ./ (c +  tau).^2;

% Convolve slew rate waveform with impulse response.
s = diff(g, 1, 2)/dt;     % T/m/s
p = zeros(3, n);
for ch = 1:3
    tmp = arg.wt(ch) * 100 * conv(s(ch,:), f);  % percent of stimulation threshold
    p(ch,:) = tmp(1:n);
end

% return total (gradient-combined) PNS waveform
pt = sqrt(sum(p.^2,1));    

% plot
if arg.plt
    tt = dt*1e3*(1:n);  % ms

    subplot(221); plot(tt, 1e3*g); ylabel('g (mT/m)'); grid on; xlabel('ms');
    axis([0 tt(end) 1.1*1e3*max(abs(g(:)))*[-1 1] ]);

    subplot(223); plot(tt(1:end-1), s); ylabel('slew (T/m/s)'); grid on; xlabel('ms');
    axis([0 tt(end) 1.1*max(abs(s(:)))*[-1 1] ]);

    subplot(224); plot(tt, p, '--'); ylabel('PNS (% of threshold)'); grid on; xlabel('ms');
    hold on; plot(tt, pt); 
    legend('x', 'y', 'z', 'total');
    axis([0 tt(end) 1.1*max(pt)*[-1 1] ]);
end

return


function p_this = sub_test
    % Validate against toppe.pns()

    dt = 4e-6;   % gradient raster time, sec
 
    % XRM coil (MR750)
    c = 360e-6;           
    rheobase = 23.4;
    alpha = 0.333;
    Smin = rheobase/alpha;

    % gradient waveform to test
    dur = 1e-3;   % sec
    ramp = linspace(0, 1, 100);
    g = [ramp ones(1,200), fliplr(ramp)];
    g = g(1:end-1);
    g = repmat([g -g], [1 5]);
    g = [g; 0.5*g; 0.5*g];
    g = 4 * g * 1e-2;    % T/m

    % calculate pns
    [p_this, ~] = pge2.pns(Smin, c, g, dt, 'plt', true);

    % ground truth
    % to compare, set wt = [1 1 1] above
    %figure;
    %p_true = toppe.pns(g, 'xrm', 'plt', true);
    %err = norm(p_true-p_this)/norm(p_true);
    %fprintf('norm(p_true-p_this)/norm(p_true)  = %.3f\n', err);
return
