function alpha = getrotationangle(b1, b2, AX)
% function R = getrotationangle(b1, b2, AX)
%
% Returns the angle alpha that rotates the gradients in b2
% to those in b1 about axis AX, if such an angle exists.
%
% Inputs
%  b1       Pulseq block containing arbitrary gradients
%  b2       Pulseq block containing arbitrary gradients
%  AX       'x', 'y', or 'z'. rotation axis
%
% Output
%  alpha   [1]  radians. Return empty if none exists.
%               

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

% zero out on-axis component
switch AX
    case 'x'
       G1(:,1) = 0;
       G2(:,1) = 0;
       u = [1 0 0];
    case 'y'
       G1(:,2) = 0;
       G2(:,2) = 0;
       u = [0 1 0];
    case 'z'
       G1(:,3) = 0;
       G2(:,3) = 0;
       u = [0 0 1];
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
u = C(I(1),:);   

% Rotation angle alpha
D = dot(G1, G2, 2)./[vecnorm(G1,2,2).*vecnorm(G2,2,2)];
alpha = acos(mode(D));  % radians. This will probably fail if waveform contains mostly zeros.

% check that rotating b2 indeed matches the gradients in b1
%if norm(G1' - R*G2', "fro")/norm(G1, "fro") > 1e-3
%    keyboard
%end

R = angleaxis2rotmat(alpha, u);
assert(norm(G1' - R*G2', "fro")/norm(G1, "fro") < 1e-3, 'Rotation detection failed');
