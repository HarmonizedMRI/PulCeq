function [wout, Ikeep] = robustinterp1(ttin, win, ttout)
% robustinterp1 - 1D interpolation
%
% function [wout, Ikeep] = robustinterp1(ttin, win, ttout)
%
% Inputs
%   ttin    [m]   
%   win     [m]
%   ttout   [n]   
%
% Output
%   wout    size(ttout)
%   Ikeep   [m]           non-NaN indices to keep in ttin/win arrays

[nr, nc] = size(ttout);  

% make times unique
ttin = ttin(:) + 1e-12*(0:numel(ttin)-1)';
ttout = ttout(:) + 1e-12*(0:numel(ttout)-1)';

% interpolate and remove NaN
w = interp1(ttin, win(:), ttout);

Ikeep = ~isnan(w);
wout = w(Ikeep);

% return as either column or row vector to match ttout
if nr > nc
    wout = wout(:);
else
    wout = wout(:).';
end


