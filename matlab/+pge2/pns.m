function [pt, p] = pns(Smin, c, g, dt, plt)
% function [pt, p] = pns(Smin, c, g, dt, [plt=false])
%
% Calculate PNS following IEC 60601-2-33:2022.
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
%   plt    true/false   plot 
%
% Output
%   pt      [1 n]  channel-combined PNS waveform (% of stimulation threshold)
%   p       [3 n]  PNS waveform on each gradient channel
%
% To test

if ischar(Smin)
    if strcmp(Smin, 'test')
        pt = sub_test();
        return;
    end
end

if nargin < 5
    plt = false;
end

assert(size(g,1) < size(g,2), 'g must be size 3xn');
assert(size(g,1) == 3, 'g must be size 3xn');

n = size(g,2);

% Contribution of a slew impulse at time 0 to PNS at time tau
% IEC 60601-2-33:2022 Eq. AA.21
tau = dt:dt:(10*c);  % no need to make much longer than chronaxie
f = dt/Smin * c ./ (c +  tau).^2;

% Convolve slew rate waveform with impulse response.
% x/y/z channel weightings from IEC 60601-2-33:2022 section (12) 
s = diff(g, 1, 2)/dt;     % T/m/s
p = zeros(3, n);
wt = [0.8 1.0 0.7];  % channel weights
for ch = 1:3
    tmp = wt(ch) * 100 * conv(s(ch,:), f);  % percent of stimulation threshold
    p(ch,:) = tmp(1:n);
end

% return total (gradient-combined) PNS waveform
pt = sqrt(sum(p.^2,1));    

% plot
if plt
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
    Smin = 20.0;

    % gradient waveform to test
    dur = 1e-3;   % sec
    n = round(dur/dt);
    g = zeros(3, n);
    t = dt*[1:n] - dt;
    g(1,:) = 0.02 * sin(2 * pi * 100 * t);
    g(2,:) = 0.015 * sin(2 * pi * 100 * t);
    g(3,:) = 0.018 * sin(2 * pi * 100 * t);
    %ramp = linspace(0, 1, 100);
    %g = [ramp ones(1,200), fliplr(ramp)];
    %g = g(1:end-1);
    %g = repmat([g -g], [1 round(n/length(g)/2)]);
    %g = [g; 0.85*g; 0.75*g];
    %g = 4 * g * 1e-2;    % T/m

    % calculate pns
    tic; p_this = pns(Smin, c, g, dt, true); ct = toc;
    fprintf('Compute time: %.4f s\n', ct);

    p_this(1:20)

    % ground truth
    %figure;
    %tic; p_true = toppe.pns(g, 'xrm', 'plt', true); ct = toc;
    %fprintf('toppe.pns() compute time: %.4f s\n', ct);
    %assert(norm(p_true-p_this)/norm(p_true) < 0.05, ...
    %    'test failed');

    return
