function d = read_cb_export(p, file)
% Validate a completed CB export before constructing its Phase-2 model.
d = jsondecode(fileread(file));
assert(strcmp(d.format, 'smiqp-miqcr-binary-v1'), 'Unknown CB export format.');
n = numel(p.c);
assert(d.n == n && d.nb_int == n && d.m == size(p.G, 1) && ...
    d.p == size(p.A, 1), 'CB export dimensions do not match the input.');
assert(all([d.mq, d.pq] >= 0) && all(mod([d.mq, d.pq], 1) == 0));

% Check the dimensions and finiteness of every coefficient used by cb_model.
fields = {'q', 'new_q', 'y_mask', 'c', 'lower_bounds', 'upper_bounds', ...
    'a', 'b', 'd', 'e', 'aq', 'bq', 'dq', 'eq', 'constant', 'new_lambda_min'};
sizes = [n*n, n*n, n*n, n, n, n, d.m*n, d.m, d.p*n, d.p, ...
    d.mq*(n+1)^2, d.mq, d.pq*(n+1)^2, d.pq, 1, 1];
for k = 1:numel(fields)
    v = d.(fields{k});
    assert(isnumeric(v) && isreal(v) && numel(v) == sizes(k) && ...
        all(isfinite(v(:))), 'Invalid CB export field: %s.', fields{k});
end
assert(all(d.lower_bounds == 0) && all(d.upper_bounds == 1));
assert(d.constant == 0, 'Unexpected constant in CB export.');
Q = reshape(d.new_q, n, n)';
assert(norm(Q-Q', 'fro') <= 1e-10*max(1, norm(Q, 'fro')), ...
    'The exported quadratic matrix is not symmetric.');
mask = reshape(d.y_mask, n, n)';
assert(all(ismember(mask(:), [0, 1])) && isequal(mask, mask'), ...
    'Invalid CB product mask.');

% The writer passes the original constant separately through cb_model.
assert(norm(reshape(d.q, n, n)'-p.Q, 'fro') < 1e-10);
assert(norm(reshape(d.d, n, d.p)'-p.A, 'fro') < 1e-10);
assert(norm(reshape(d.a, n, d.m)'-p.G, 'fro') < 1e-10);
assert(norm(d.c(:)-p.c(:))+norm(d.e(:)-p.b(:))+ ...
    norm(d.b(:)-p.h(:)) < 1e-10, 'CB export does not match the input.');
end
