function [tt, wav] = trap2arb(g)
% function [tt, wav] = g2arb(g)
%
% Convert trapezoid to arbitrary gradient samples at corner points.
% First sample time is g.delay.
%
% Inputs
%   g      struct      Pulseq gradient event
%
% Ouputs
%   tt    [n 1]        sample time points ('corner' points), with tt(1) = g.delay
%   wav   [n 1]        gradient samples at tt

if g.flatTime > 0
    wav = [0 1 1 0]*g.amplitude;       
    tt = g.delay + [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime];
else
    wav = [0 1 0]*g.amplitude;       
    tt =  g.delay + [0 g.riseTime g.riseTime+g.fallTime];
end

