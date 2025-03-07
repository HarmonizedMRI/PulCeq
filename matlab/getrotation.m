function [R,scale] = getrotation(b1, b2)
% function [R,scale] = getrotation(b1, b2)
%
% Determine if a rotation matrix R exists that takes (arbitrary) gradients from
% Pulseq block b2 to those in b1.
%
% For trapezoids, just return identity since rotation just
% results in a trapezoid scaling on each axis which seq2ceq already detects.
%
% Inputs
%  b1   Pulseq block containing arbitrary gradients
%  b2   Pulseq block containing arbitrary gradients
%
% Output
%  R       [3 3] or []    Rotation matrix taking gradients in b2 to those in b1.
%                         If none is found return empty matrix.
%  scale   [1] or []      Scaling factor from gradients in b2 to b1. 
%                         If no rotation found, return empty.

% Get number of samples N
N = 0;
N2 = 0;
tt  = [];   % sample times
tt2 = [];
for ax = {'gx','gy','gz'}
    g = b1.(ax{1});

    % If trapezoid, just return identity
    if isfield(g, 'type')
        if strcmp(g.type, 'trap')
            scale = 1.0;
            R = eye(3);
            return;
        end
    end

    g2 = b2.(ax{1});

    if isfield(g, 'waveform')
        N = max(N, length(g.waveform));
        tt = g.tt;
    end

    if isfield(g2, 'waveform')
        N2 = max(N2, length(g2.waveform));
        tt2 = g2.tt;
    end
end

% Insist the sample times and number of waveform samples are the same 
if ~all(abs(tt-tt2) < eps) | N ~= N2
    R = [];
    scale = [];
    return;
end

% Construct 3D vectors 
nd = 3;
G1 = zeros(N, nd);
G2 = zeros(N, nd);
ax = {'gx','gy','gz'};
for ii = 1:nd
    g1 = b1.(ax{ii});
    if isfield(g1, 'waveform')
        G1(:,ii) = g1.waveform;
    end
    g2 = b2.(ax{ii});
    if isfield(g2, 'waveform')
        G2(:,ii) = g2.waveform;
    end
end

% Solve for rotation using Procruste's analysis:
% i.e. for given set of points A and B, where B is a unitary transformation
% of A (no scaling/translation),
% solve for unitary matrix Q* = argmin_Q ||QA - B||^2
% --> closed form solution: Q* = UV' where BA' = USV'
% let B = G1', A = G2':
[U,S,V] = svd(G1' * G2);
R = U*V';

% If R is improper (negative determinant), fix it by negating last RSV
if det(R) < 0
    V(:, end) = -V(:, end);
    R = U * V';
end

% Determine scale factor
% after finding Q* which is a scale invariant solution,
% solve for scaling factor alpha* = argmin_q ||alphaQ*A - B||^2
% --> closed form solution: alpha* = trace(S)/trace(A'A)
scale = trace(S)/norm(G2','fro')^2;

% check that rotating b2 indeed matches the gradients in b1, otherwise
% assume b1 is not a rotation of b2
if norm(G1' - scale*R*G2', "fro")/norm(G1, "fro") > 1e-3
    R = [];
    scale = [];
end

