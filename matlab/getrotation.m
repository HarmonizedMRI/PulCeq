function R = getrotation(b1, b2)
% function R = getrotation(b1, b2)
%
% Determine if a rotation matrix R exists that takes (arbitrary) gradients from
% block b1 to those in b2.
%
% Note that if a trapezoid is rotated, this function isn't necessary since
% rotation just results in a trapezoid scaling which seq2ceq already detects.
%
% Inputs
%  b1   Pulseq block containing arbitrary gradients
%  b2   Pulseq block containing arbitrary gradients
%
% Output
%  R    [3 3] or []    Rotation matrix taking gradients in b1 to those in b2.

% Get number of samples.
% If not the same, return []
N = 0;
N2 = 0;
tt  = [];   % sample times
tt2 = [];
for ax = {'gx','gy','gz'}
    g = b1.(ax{1});
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

if ~all(tt == tt2) | N ~= N2
    R = [];
    return;
end

% Construct vectors 
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

% Get rotation axis (cross product)
C = zeros(N,nd);
for n = 1:N
    C(n,:) = cross(G1(n,:), G2(n,:));
end

%if all(diff(C(n)) < eps) & norm(
%    'success'
%end

