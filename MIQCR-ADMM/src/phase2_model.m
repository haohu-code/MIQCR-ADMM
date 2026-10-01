function model = phase2_model(p, D, M, z, full_rlt, support)
% Form z+[1;x]'D[1;x]+<M,Y>, with binary x and continuous product variables.

% Select one variable for each retained off-diagonal product x_i*x_j.
n = numel(p.c);
if full_rlt
    [i, j] = find(triu(ones(n), 1));
elseif nargin == 6
    [i, j] = find(triu(support, 1));
else
    [i, j] = find(triu(M(2:end, 2:end) ~= 0, 1));
end
k = numel(i);
I = sparse((1:k)', i, 1, k, n);
J = sparse((1:k)', j, 1, k, n);

% Keep D's quadratic block; express lifted M terms linearly in x and products.
model.Q = blkdiag(sparse(D(2:end, 2:end)), sparse(k, k));
model.obj = [2*D(1, 2:end)'+2*M(1, 2:end)'+diag(M(2:end, 2:end)); ...
    2*M(sub2ind([n+1, n+1], i+1, j+1))];
model.objcon = z+D(1, 1)+M(1, 1);

% Original constraints and McCormick links: y<=x_i, y<=x_j, y>=x_i+x_j-1.
model.A = [p.A, sparse(size(p.A, 1), k); p.G, sparse(size(p.G, 1), k); ...
    -I, speye(k); -J, speye(k); I+J, -speye(k)];
model.rhs = [p.b; p.h; zeros(2*k, 1); ones(k, 1)];
model.sense = [repmat('<', size(p.A, 1), 1); ...
    repmat('=', size(p.G, 1), 1); repmat('<', 3*k, 1)];
model.lb = zeros(n+k, 1);
model.ub = ones(n+k, 1);
model.vtype = [repmat('B', n, 1); repmat('C', k, 1)];
model.modelsense = 'min';

if full_rlt
    % Optional lifted constraints; normalization/arrow are built into the map.
    % vec(Y)=e00+T*[x;y].
    d = n+1;
    T = product_map(n, i, j);
    H = [p.b, -p.A; ones(n, 1), -speye(n)];
    R = kron(speye(d), H)*T;
    c = [H(:, 1); zeros(n*size(H, 1), 1)];
    rows = size(H, 1);
    spoke = rows+1:d*rows;
    cap = repmat((1:rows)', n, 1);
    model.A = [model.A; -R(spoke, :); R(spoke, :)-R(cap, :)];
    model.rhs = [model.rhs; c(spoke); c(cap)-c(spoke)];
    model.sense = [model.sense; repmat('<', 2*n*rows, 1)];
    F = [p.h, -p.G];
    model.A = [model.A; kron(speye(d), F)*T];
    model.rhs = [model.rhs; -F(:, 1); zeros(n*size(F, 1), 1)];
    model.sense = [model.sense; repmat('=', d*size(F, 1), 1)];
end
model.rhs = full(model.rhs);
model.obj = full(model.obj);
end
