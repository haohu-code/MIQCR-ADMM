function r = equality_admm(p, tol, maxit, limit, maxchanges)
% Solve the equality-face DNN relaxation by PSD and polyhedral projections.
% Return the finite iterate and multiplier estimates for later reformulation.
if nargin < 4
    limit = Inf;
end
if nargin < 5
    maxchanges = 8;  % Zero gives a fixed penalty throughout.
end
assert(isempty(p.b), 'Use row-star ADMM for inequality constraints.');
clock = tic;
status = 'ITERATION_LIMIT';

% Restrict the lifted matrix to the nullspace of the equality system.
n = numel(p.c);
F = [p.h, -p.G];
V = null(full(F));
P = V*V';
C_original = [p.constant, p.c'/2; p.c/2, p.Q];
C = C_original;
qap = strcmp(p.family, 'qap');
box = isfield(p, 'box_gangster') && p.box_gangster;
if ~qap
    C = P*C*P;
end
scale = n/max(norm(C, 'fro'), eps);
Y = p.Y0;
R = V'*Y*V;
Z = Y-V*R*V';

% Adapt the penalty only finitely many times, then keep it fixed.
if qap
    beta = sqrt(n)/3;
    omega = 1;
    stable_needed = 5;
else
    beta = n/3;
    omega = 1.618;
    stable_needed = 2;
end
stable = 0;
changes = 0;
last_change = 0;
for iter = 1:maxit
    % PSD projection in the reduced space, followed by the polyhedral step.
    old = Y;
    R = psd(V'*(Y+Z/beta)*V);
    W = V*R*V'-(scale*C+Z)/beta;
    W = (W+W')/2;
    if qap
        % The arrow equations follow at consensus from assignment structure.
        Y = min(1, max(0, W));
        Y(p.gangster) = 0;
        Y(1, 1) = 1;
    else
        Y = max(0, W);
        if box
            Y = min(1, Y);
        end
        Y(1, 1) = 1;
        for i = 2:n+1
            % Project the diagonal and its two symmetric arrow entries jointly.
            Y(i, i) = max(0, (2*W(1, i)+W(i, i))/3);
            if box
                Y(i, i) = min(1, Y(i, i));
            end
            Y(1, i) = Y(i, i);
            Y(i, 1) = Y(i, i);
        end
    end
    if box
        Y(p.gangster) = 0;
    end

    % Multiplier update and relative primal/dual residuals.
    K = beta*(Y-W)/scale;
    Z = Z+omega*beta*(Y-V*R*V');
    Z = (Z+Z')/2;
    a = norm(Y-V*R*V', 'fro');
    b = beta*norm(V'*(Y-old)*V, 'fro');
    u = max(norm(Y, 'fro'), norm(R, 'fro'));
    v = norm(V'*Z*V, 'fro');
    residual = max(a/max(1, u), b/max(1, v));
    if residual < tol
        stable = stable+1;
    else
        stable = 0;
    end
    if stable >= stable_needed
        status = 'CONVERGED';
        break;
    end
    if toc(clock) >= limit
        status = 'TIME_LIMIT';
        break;
    end
    if changes < maxchanges && mod(iter, 25) == 0 && iter < maxit
        next = beta;
        if a/max(u, eps) > 5*b/max(v, eps)
            next = 2*beta;
        elseif b/max(v, eps) > 5*a/max(u, eps)
            next = beta/2;
        end
        if next ~= beta
            % Z is an unscaled multiplier and must not be rescaled.
            beta = next;
            changes = changes+1;
            last_change = iter;
        end
    end
end

% Undo objective projection so C+Z uses the original objective convention.
r = struct('Y', Y, 'Z', Z/scale+C-C_original, 'Znative', Z/scale, ...
    'Cused', C, 'V', V, 'C', C_original, 'iter', iter, ...
    'residual', residual, 'primal_residual', a/max(1, u), ...
    'dual_residual', b/max(1, v), 'status', status, ...
    'seconds', toc(clock), 'penalty', beta, 'penalty_updates', changes);
end

function X = psd(X)
% Orthogonal projection onto the positive semidefinite cone.
[U, e] = eig((X+X')/2, 'vector');
X = (U.*max(e, 0)')*U';
X = (X+X')/2;
end
