function check = validate_incumbent(p, model, answer, pairs, error_bound)
% Check the original binary point and its exactly reconstructed products.
% Small deviations in the solver vector are recorded, not confused with a
% failure of the reformulation identity at the reconstructed binary point.
n = numel(p.c);
raw = answer.x(:);
assert(numel(raw) == numel(model.obj) && all(isfinite(raw)));
x = round(raw(1:n));
integrality = max(abs(raw(1:n)-x));
assert(integrality < 1e-5 && all(ismember(x, [0, 1])));
assert(norm(p.G*x-p.h, Inf) < 1e-7, 'Original equalities violated.');
v = [x; x(pairs(:,1)).*x(pairs(:,2))];
residual = model.A*v-model.rhs;
assert(all(abs(residual(model.sense == '=')) < 1e-7) && ...
    all(residual(model.sense == '<') < 1e-7) && ...
    all(residual(model.sense == '>') > -1e-7), ...
    'Reconstructed point violates model constraints.');
assert(all(v >= model.lb-1e-7) && all(v <= model.ub+1e-7));
original = x'*p.Q*x+p.c'*x+p.constant;
reconstructed = v'*model.Q*v+model.obj'*v+model.objcon;
raw_value = raw'*model.Q*raw+model.obj'*raw+model.objcon;
tolerance = 1e-5+1e-9*abs(original);
assert(abs(reconstructed-original) <= error_bound+tolerance, ...
    'Binary reformulation objective mismatch.');
assert(isfinite(answer.objval) && abs(raw_value-answer.objval) <= tolerance, ...
    'Reported solver objective disagrees with its returned vector.');
check = struct('x', x, 'objective', original, ...
    'binary_identity_error', reconstructed-original, ...
    'raw_solver_error', raw_value-answer.objval, ...
    'rounding_objective_change', reconstructed-raw_value, ...
    'integrality_deviation', integrality);
end
