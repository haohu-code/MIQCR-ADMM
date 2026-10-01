function [model, error_bound] = cb_model(p, d)
% Translate the SMIQP export into a convex MIQP with one copy of each product.

% Separate the convex quadratic block from the retained product coefficients.
n = d.n;
Q = reshape(d.new_q, n, n)'/2-d.new_lambda_min*eye(n);
K = reshape(d.q, n, n)'-Q;
mask = reshape(d.y_mask, n, n)' ~= 0;
% Since products lie in [0,1], this bounds the error from omitted coefficients.
error_bound = sum(abs(K(~mask)));
K(~mask) = 0;
D = [d.constant, d.c'/2; d.c/2, Q];
M = blkdiag(0, K);

% Transfer a diagonal shift between blocks without changing binary values.
shift = max(0, -min(eig((Q+Q')/2)))+1e-9;
D(2:end, 2:end) = Q+shift*eye(n);
M(2:end, 2:end) = K-shift*eye(n);
model = phase2_model(p, D, M, 0, false, mask);
if isfield(p, 'constant')
    model.objcon = model.objcon+p.constant;
end
[i, j] = find(triu(mask, 1));
T = product_map(n, i, j);

% Keep the linearizations of every quadratic constraint in the export.
for r = 1:d.mq+d.pq
    if r <= d.mq
        v = d.aq((r-1)*(n+1)^2+(1:(n+1)^2));
        rhs = d.bq(r);
        sense = '=';
    else
        k = r-d.mq;
        v = d.dq((k-1)*(n+1)^2+(1:(n+1)^2));
        rhs = d.eq(k);
        sense = '<';
    end
    A = reshape(v, n+1, n+1)';
    assert(~any(triu((A(2:end, 2:end)+A(2:end, 2:end)') ~= 0, 1) & ~mask, 'all'), ...
        'The export omits a product required by a quadratic constraint.');
    model.A(end+1, :) = A(:)'*T;
    model.rhs(end+1, 1) = rhs-A(1, 1);
    model.sense(end+1, 1) = sense;
end
end
