function [model, d] = continuous_bound_model(p, r)
% Construct an exact convex MIQP from finite ADMM multiplier estimates.
% Optimize the continuous bound with the centered arrow and product terms fixed.
% Here p.G*x = p.h is the equality system denoted Ax = b in the paper;
% p.A*x <= p.b stores any additional linear inequalities.

% Extract multipliers from the stationarity estimate C+Z.
n = numel(p.c);
T = r.C+r.Z;
T = (T+T')/2;
if strcmp(p.family, 'qap') || (isfield(p, 'box_gangster') && p.box_gangster)
    % Retain signed box terms; gangster terms vanish at feasible assignments.
    assert(isequal(size(p.gangster), size(T)) && isequal(p.gangster, p.gangster'));
    G = zeros(n+1);
    G(p.gangster) = T(p.gangster);
    T(p.gangster) = 0;
    T(1, 1) = 0;
    N = max(T, 0);
    U = max(-T, 0);
    M = N-U;
    bar = zeros(n, 1);
    S = r.C-M-G;
else
    % Frobenius projection onto arrow terms plus a nonnegative lower block.
    N = max(0, T(2:end, 2:end));
    lambda = zeros(n, 1);
    for i = 1:n
        q = T(i+1, i+1)+2*T(1, i+1);
        if q >= 0
            lambda(i) = 2*T(1, i+1);
            N(i, i) = q;
        else
            lambda(i) = 2/3*(T(1, i+1)-T(i+1, i+1));
            N(i, i) = 0;
        end
    end
    bar = lambda-mean(lambda);
    M = blkdiag(0, N);
    S = r.C-[0, bar'/2; bar/2, -diag(bar)]-M;
end

% Choose the smallest shift giving nonnegative curvature on null(A).
A = full(p.G);
W = null(A);
P = W*W';
K0 = S(2:end, 2:end);
H = W'*K0*W;
mu = 0;
if ~isempty(H)
    mu = -min(eig((H+H')/2));
end
K = K0+mu*eye(n);

% Complete the quadratic block using terms that vanish on Ax=b.
if isempty(A)
    Gamma = zeros(0, n);
else
    Gamma = -(A*A')\(A*K*(eye(n)+P));
end
Q = P*K*P;
Q = (Q+Q')/2;
c = 2*S(2:end, 1)-mu*ones(n, 1)-Gamma'*p.h;
D = [S(1, 1), c'/2; c/2, Q];
model = phase2_model(p, D, M, 0, false);

% Add epsilon*sum(x_i^2-x_i) for a numerical PSD margin, preserving binary values.
e = min(eig(Q));
guard = max(0, 1e-8*max(1, norm(Q, Inf))-e);
model.Q(1:n, 1:n) = model.Q(1:n, 1:n)+guard*speye(n);
model.obj(1:n) = model.obj(1:n)-guard;
recovery = norm(Q-K-(A'*Gamma+Gamma'*A)/2, 'fro')/max(1, norm(Q, 'fro'));
assert(recovery < 1e-7 && e+guard >= 0, 'Convex completion failed.');
d = struct('mu', mu, 'guard', guard, 'recovery', recovery, 'mineig', e+guard, ...
    'M', M, 'lambda', bar+mu, 'Gamma', Gamma);
end
