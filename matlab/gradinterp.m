function [wav, tt] = gradinterp(g, sysGE, varargin)
% function wav = gradinterp(g, sysGE, varargin)
%
% Interpolate gradient waveforms and convert to Gauss/cm
%
% We assume that there are 3 types of gradients:
% 1: arbitrary gradient specified on regular raster time. Assumed to start and end at zero ('first' and 'last' points). 
% 2: extended trapezoid, specified on "corners" of shape. Assumed to start at zero if and only if delay > 0.
% 3: trapezoid, starting and ending at 0
%
% To run test function, do:
% >> gradinterp('test');
%
% Inputs
%    g         struct    Pulseq gradient event
%    sysGE     struct    GE hardware settings, see toppe.systemspecs()
% 
% Output
%    wav   [n 1]     Gradient waveform interpolated to 4us raster time, Gauss/cm

if isstr(g)
    sub_test();
    return;
end

% parse input options
arg.seqGradRasterTime = 10e-6;
arg.outRasterTime = 4e-6;
arg = vararg_pair(arg, varargin);

gamma = sysGE.gamma;   % Hz/T
raster = sysGE.raster*1e-6;   % sec

if isempty(g)
    wav = [];
else
    if strcmp(g.type, 'grad')
        % Arbitrary gradient
        % restore shape: if we had a
        % trapezoid converted to shape we have to find
        % the "corners" and we can eliminate internal
        % samples on the straight segments
        % but first we have to restore samples on the
        % edges of the gradient raster intervals
        % for that we need the first sample

        % check if we have an extended trapezoid, or an arbitrary gradient on a regular raster
        tt_rast=g.tt/arg.seqGradRasterTime+0.5;
        if all(abs(tt_rast-(1:length(tt_rast))')<1e-6)  % samples assumed to be on center of raster intervals
            % arbitrary gradient on a regular raster
            %g.first = 0;
            %g.last = 0;
            areaIn = sum(g.waveform)*arg.seqGradRasterTime;
            wavtmp = [g.first g.waveform(:)' g.last];
            tttmp = g.delay + [0 g.tt(:)' g.tt(end)+arg.seqGradRasterTime/2];
        else
            % extended trapezoid: shape specified on "corners" of waveform
            areaIn = sum( (g.waveform(1:(end-1)) + g.waveform(2:end))/2 .* diff(g.tt) );
            wavtmp = g.waveform(:)';
            tttmp = g.delay + g.tt(:)';
        end
    else
        % Convert trapezoid to arbitrary gradient
        areaIn = [g.riseTime/2 + g.flatTime + g.fallTime/2] * g.amplitude;
        [tttmp, wavtmp] = trap2arb(g);
    end

    % Add delay
    if g.delay > 0
        wavtmp = [0 wavtmp];
        tttmp = [0 tttmp];
    end

    % interpolate (includes the delay)
    tt = raster/2 : raster : tttmp(end);
    tmp = interp1(tttmp, wavtmp, tt);

    areaOut = sum(tmp) * raster;

    if any(isnan(tmp))
        msg = sprintf('NaN in gradient trapezoid waveform after interpolation (parent block %d)', p);
        error(msg);
    end

    % If areas don't match to better than 0.01%, throw warning
    if abs(areaIn) > 1e-6 
        if abs(areaIn-areaOut)/abs(areaIn) > 1e-4
            msg = sprintf('Gradient area not preserved after interpolating to GE raster time (in: %.3f 1/m, out: %.3f). Did you wrap all gradient events in trap4ge()?', areaIn, areaOut);
            warning(msg);
        end
    end

    % convert to Gauss/cm
    wav = tmp(:).' / gamma * 100; % Gauss/cm
end

return

function sub_test()

sys = mr.opts();   % default settings
fov = 24e-2;       % m
deltak = 1/fov;

% define two trapezoids: before and after wrapping in trap4ge
g = mr.makeTrapezoid('x', sys, 'Area', 1*deltak);
g2 = trap4ge(g, 20e-6, sys);  % ensure that sample points are on 20us boundary
g.delay = 40e-6;              % on 20us boundary
g2.delay = 40e-6;             % on 20us boundary

subplot(131); 
title('Trapezoid interpolated to 4us, before wrapping in trap4ge');
sub_test_doone(g);

subplot(132);
title('Trapezoid interpolated to 4us, after wrapping in trap4ge');
sub_test_doone(g2);

subplot(133);
title('Plotted together');
sub_test_doone(g);
sub_test_doone(g2);

return

function sub_test_doone(g)

sysGE = toppe.systemspecs();   % default settings
raster = sysGE.raster*1e-6;   % sec

wav = gradinterp(g, sysGE);
tt = (1:length(wav))*raster - raster/2;  % us

% plot in units of Gauss/cm
g.amplitude = g.amplitude / sysGE.gamma * 100;   % G/cm
[ttin, wavin] = trap2arb(g);

hold on;
plot(ttin*1e3,wavin,'b-x'); 
plot(tt*1e3,wav,'ro'); 
xlabel('time (ms)');
ylabel('Gauss/cm');

areaIn = (g.riseTime/2 + g.flatTime + g.fallTime/2)*g.amplitude;
areaOut = sum(wav) * raster;

legend(sprintf('In: area = %.4f G/cm*us', areaIn*1e6) , ...
       sprintf('Out: area = %.4f G/cm*us', areaOut*1e6));

fprintf('Areas = %.5f (in), %.5f (out) G/cm*us\n', areaIn*1e6, areaOut*1e6);

return
