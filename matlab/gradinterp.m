function wav = gradinterp(g, sysGE, varargin)
% Interpolate gradient waveforms and convert to Gauss/cm
% We assume that there are 3 types of gradients:
% gradType = 1: arbitrary gradient specified on regular raster time. Assumed to start and end at zero ('first' and 'last' points). 
% gradType = 2: extended trapezoid, specified on "corners" of shape. Assumed to start at zero if and only if delay > 0.
% gradType = 3: trapezoid, starting and ending at 0
%
% Inputs
%    g         struct    Pulseq gradient event
%    sysGE     struct    GE hardware settings, see toppe.systemspecs()
% 
% Output
%    wav   [n 1]     Gradient waveform interpolated to 4us raster time, Gauss/cm

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
            gradType = 1;
            areaIn = sum(g.waveform)*arg.seqGradRasterTime;
            g.first = 0;
            g.last = 0;
            wavtmp = [g.first g.waveform(:)' g.last];
            tttmp = g.delay + [0 g.tt(:)' g.tt(end)+arg.seqGradRasterTime/2];
        else
            % extended trapezoid: shape specified on "corners" of waveform
            gradType = 2;
            areaIn = sum( (g.waveform(1:(end-1)) + g.waveform(2:end))/2 .* diff(g.tt) );
            wavtmp = g.waveform(:)';
            tttmp = g.delay + g.tt(:)';
        end
    else
        % Convert trapezoid to arbitrary gradient
        gradType = 3;
        areaIn = [g.riseTime/2 + g.flatTime + g.fallTime/2] * g.amplitude;
        if g.flatTime > 0
            wavtmp = [0 1 1 0]*g.amplitude;       
            tttmp = g.delay + [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime ];
        else
            wavtmp = [0 1 0]*g.amplitude;       
            tttmp = g.delay + [0 g.riseTime g.riseTime+g.fallTime ];
        end
    end

    if g.delay > 0
        wavtmp = [0 wavtmp];
        tttmp = [0 tttmp];
    end

    tge = raster/2 : raster : tttmp(end);
    tmp = interp1(tttmp, wavtmp, tge);

    if any(isnan(tmp))
        msg = sprintf('NaN in gradient trapezoid waveform after interpolation (parent block %d)', p);
        error(msg);
    end

    % convert to Gauss/cm
    wav = tmp(:).' / gamma * 100; % Gauss/cm
end
