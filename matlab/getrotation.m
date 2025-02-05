function R = getrotation(b1, b2)
% function R = getrotation(b1, b2)
%
% Determine if a rotation matrix R exists that takes (arbitrary) gradients from
% Pulseq block b2 to those in b1.
%
% For trapezoids, this function isn't necessary since rotation just
% results in a trapezoid scaling on each axis which seq2ceq already detects.
%
% Inputs
%  b1   Pulseq block containing arbitrary gradients
%  b2   Pulseq block containing arbitrary gradients
%
% Output
%  R    [3 3] or []    Rotation matrix taking gradients in b2 to those in b1.
%                      If none is found return empty matrix.

% Get number of samples N
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

% Insist the sample times and number of waveform samples are the same 
if ~all(abs(tt-tt2) < eps) | N ~= N2
    R = [];
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

% Get rotation axis for each sample time
C = cross(G2, G1, 2);

% If all cross products are zero, return identity
% This also works if C = [] (i.e., block has no gradients)
if rank(C) < 100*eps
    R = eye(3);
    return;
end

% If rotation axis not the same for all samples, return []
if abs(rank(C)-1) > 100*eps
    R = [];  
    return;
end

% Axis of rotation u
A = vecnorm(C');
I = find(A == max(A));  % avoid samples with zero gradient amplitude
u = C(I,:);   

% Rotation angle alpha
D = dot(G1, G2, 2)./[vecnorm(G1,2,2).*vecnorm(G2,2,2)];
alpha = acos(mode(D));  % radians. This will probably fail if waveform contains mostly zeros.

R = angleaxis2rotmat(alpha, u);

% check that rotating b2 indeed matches the gradients in b1
assert(norm(G1' - R*G2', "fro")/norm(G1, "fro") < 1e-3, 'Rotation detection failed');
