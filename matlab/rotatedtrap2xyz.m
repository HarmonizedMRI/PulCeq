function [gx, gy, gz] = rotatedtrap2xyz(g)
% function [gx, gy, gz] = rotatedtrap2xyz(g)
%
% Input
%  g    length-1/length-2/length-3 cell array of gradient trapezoids, as returned by mr.rotate()
%
% Output
%  gx/gy/gz  x/y/z gradient events. An empty/zero-amplitude gradient is given amplitude = eps
%            to preserve it as a trapezoid in the .seq file.

% Initialize with zero amplitude
for ax = {'x', 'y', 'z'}
    eval(sprintf('g%s = g{1};', ax{1}));
    eval(sprintf("g%s.channel = '%s';", ax{1}, ax{1}));
    eval(sprintf("g%s.amplitude = eps;", ax{1}));
    eval(sprintf("g%s.area = eps;", ax{1}));
    eval(sprintf("g%s.flatArea = eps;", ax{1}));
    if strcmp(eval(sprintf('g%s.type', ax{1})), 'grad')
        eval(sprintf("g%s.tt = g{1}.tt;", ax{1}));
        eval(sprintf("g%s.waveform = 0*g{1}.waveform;", ax{1}));
    end
end

% Replace with any non-zero trapezoids
for i = 1:length(g)
    eval(sprintf('g%s = g{i};', g{i}.channel));
end
