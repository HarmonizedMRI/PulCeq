function out = comple(idx, n)
% COMPLE  Complementary indices.
%   OUT = COMPLE(IDX, N) returns all indices from 1..N except those in IDX.
%
% Example:
%   comple([2 4], 6)   ->   [1 3 5 6]
%
% Notes:
%   - IDX may be a scalar or a vector.
%   - Output is always a row vector.
%   - Values outside 1..N are ignored.

    % Ensure column vector for safety
    idx = idx(:);

    % Keep only valid indices
    idx = idx(idx >= 1 & idx <= n);

    % Complement
    out = setdiff(1:n, idx, 'stable');
end

