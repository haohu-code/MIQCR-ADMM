function row = run_experiment(task, total_limit, p1_limit, folder)
% Run one instance with MIQCR-FR (using ADMM), MIQCR-CB, or direct Gurobi.
% Task = 3*(manifest row - 1) + method index, in the order ADMM, CB, GUROBI.
% ADMM is the saved method identifier for the MIQCR-FR implementation.
% P1 includes the SDP solve, reformulation, and starting-point checks.
% P2 is one direct mixed-integer solve and validation of the returned solution.

% 1. Load the instance and set solver parameters.
root = setup;
if nargin < 2
    total_limit = 3600;
end
if nargin < 3
    p1_limit = 900;
end
if nargin < 4
    folder = fullfile(root, 'outputs/continuous_bound');
end
if ~startsWith(folder, filesep)
    folder = fullfile(root, folder);
end
if ~isfolder(folder)
    mkdir(folder);
end

validateattributes(total_limit, {'numeric'}, {'scalar','finite','positive'});
validateattributes(p1_limit, {'numeric'}, {'scalar','finite','positive'});
p1_limit = min(p1_limit, total_limit);

cases = readtable(fullfile(root, 'data/manifest.csv'), ...
    'TextType', 'string', 'Delimiter', ',', ...
    'NumHeaderLines', 0, 'ReadVariableNames', true);
validateattributes(task, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', 3*height(cases)});
k = ceil(task/3);
methods = ["ADMM", "CB", "GUROBI"];
method = methods(mod(task-1, 3)+1);
s = load(fullfile(root, char(cases.File(k))));
p = s.p;
n = numel(p.c);
x0 = p.feasible;
assert(isempty(p.b), 'This runner requires equality-constrained inputs.');

file = fullfile(folder, [p.name, '_', char(method)]);
assert(~isfile([file, '.mat']) && ~isfile([file, '.csv']) && ...
    ~isfile([file, '.log']), 'Refuse to overwrite result.');
assert(all(ismember(x0, [0, 1])) && norm(p.G*x0-p.h, Inf) < 1e-7);
if isfield(s, 'incidence')
    validate_qspp_path(s, x0);
end

params = struct('OutputFlag', 1, 'Threads', 1, 'Seed', 0, ...
    'MIPGap', 1e-6, 'MIPGapAbs', 1e-6, ...
    'NonConvex', 0, 'TimeLimit', total_limit, 'LogFile', [file, '.log']);
options = params;
original = struct('Q', sparse(p.Q), 'obj', p.c, 'objcon', p.constant, ...
    'A', sparse(p.G), 'rhs', full(p.h), 'sense', repmat('=', numel(p.h), 1), ...
    'lb', zeros(n, 1), 'ub', ones(n, 1), ...
    'vtype', repmat('B', n, 1), 'start', x0);

% Initialize records so limits and failures can also be saved.
p1 = 0;
post = 0;
p2 = NaN;
bound = -Inf;
value = x0'*p.Q*x0 + p.c'*x0 + p.constant;
status = "NOT_RUN";
native = "NA";
solverstatus = "NOT_RUN";
checks = [];
pairs = zeros(0, 2);
r = [];
d = [];
model = [];
answer = [];
failure = "";
error_bound = 0;
iterations = NaN;
residual = NaN;
shift = NaN;
recovery = NaN;
mineig = NaN;
products = NaN;
nodes = 0;

% The timing marker distinguishes SDP work (1), construction (2), and P2 (3).
phase = 1;
t = tic;
try
    % 2. Phase 1: construct the convex reformulation.
    if method == "ADMM"
        r = equality_admm(p, 1e-3, 10000, p1_limit, 8);
        p1 = toc(t);
        native = string(r.status);
        iterations = r.iter;
        residual = r.residual;

        % Choose the scalar shift to optimize the continuous MIQCR bound.
        phase = 2;
        t = tic;
        [model, d] = continuous_bound_model(p, r);
        post = toc(t);
        shift = d.mu;
        recovery = d.recovery;
        mineig = d.mineig;
    elseif method == "CB"
        [d, p1] = cb_phase1(p, fullfile(folder, 'cb'), p1_limit, total_limit);
        if isempty(d)
            native = "TIME_LIMIT";
            status = "PHASE1_LIMIT";
        else
            % A validated timeout export also proceeds to Phase 2. Its full
            % elapsed time is charged to P1 before assigning the remainder.
            native = string(d.termination);
            phase = 2;
            t = tic;
            [model, error_bound] = cb_model(p, d);
            post = toc(t);
        end
    else
        model = original;
    end

    if ~isempty(model)
        % 3. Set and verify the common feasible start; charge this work to P1.
        phase = 2;
        t = tic;
        products = numel(model.obj)-n;
        if method ~= "GUROBI"
            if method == "ADMM"
                [i, j] = find(triu(d.M(2:end, 2:end) ~= 0, 1));
            else
                [i, j] = find(triu(reshape(d.y_mask, n, n)' ~= 0, 1));
            end
            pairs = [i, j];
            assert(size(pairs, 1) == products);
            model.start = [x0; x0(i).*x0(j)];
            start_answer = struct('x', model.start, ...
                'objval', model.start'*model.Q*model.start + ...
                model.obj'*model.start + model.objcon);
            validate_incumbent(p, model, start_answer, pairs, error_bound);
        end
        post = post + toc(t);

        % 4. Phase 2: submit the mixed-integer model directly to Gurobi.
        phase = 3;
        t = tic;
        remaining = total_limit-p1-post;
        if remaining <= 0
            status = "TOTAL_LIMIT";
            p2 = 0;
        else
            options = params;
            if method == "GUROBI"
                options.NonConvex = 2;
            else
                % Preserve the convex quadratic relaxation constructed in P1.
                options.PreQLinearize = 0;
            end
            % Presolve, cuts, and the root algorithm otherwise stay automatic.
            options.TimeLimit = max(0, remaining-toc(t));
            answer = gurobi(model, options);
            solverstatus = string(answer.status);
            status = solverstatus;

            % 5. Validate bounds and solutions against the original problem.
            if isfield(answer, 'objbound')
                % CB coefficient omissions require a conservative correction.
                assert(~isnan(answer.objbound) && answer.objbound < Inf);
                bound = answer.objbound-error_bound;
            end
            if isfield(answer, 'nodecount')
                nodes = answer.nodecount;
            end
            if isfield(answer, 'x')
                checks = validate_incumbent(p, model, answer, pairs, error_bound);
                if isfield(s, 'incidence')
                    validate_qspp_path(s, checks.x);
                end
                value = min(value, checks.objective);
            end
            assert(bound <= value + 1e-5 + 1e-9*abs(value), ...
                'Bound exceeds verified feasible objective.');
            p2 = toc(t); % Include post-solve validation consistently.
        end
    end
catch err
    % Retain elapsed time and diagnostic evidence when a check or solver fails.
    if phase == 1
        p1 = toc(t);
    elseif phase == 2
        post = post + toc(t);
    else
        p2 = toc(t);
    end
    status = "ERROR";
    failure = string(getReport(err, 'extended', 'hyperlinks', 'off'));
    fprintf(2, '%s\n', failure);
end

% Summarize performance, then save both the summary and the native outputs.
total = p1 + post;
if isfinite(p2)
    total = total + p2;
end
gap = 100*max(0, value-bound)/max(1, abs(value));
if status == "ERROR"
    gap = NaN;
end
threshold = max(1e-6, 1e-6*abs(value));
if p.integer_objective
    threshold = 1-1e-4;
end
proved = status ~= "ERROR" && isfinite(bound) && value-bound < threshold;
row = table(task, string(p.name), string(p.family), method, n, ...
    size(p.G, 1), n+1-size(p.G, 1), ...
    p1, post, p1+post, p2, total, value, bound, gap, proved, status, native, ...
    iterations, residual, shift, recovery, mineig, products, nodes, ...
    solverstatus, error_bound, p.integer_objective, ...
    'VariableNames', {'Task', 'Instance', 'Family', 'Method', 'N', ...
    'EqualityRank', 'ReducedOrder', 'SDPSeconds', 'Post', 'P1', 'P2', 'Total', ...
    'Objective', 'Bound', 'Gap', 'Proved', 'Status', 'Phase1Status', ...
    'Iterations', 'Residual', 'Mu', 'RecoveryResidual', 'MinEig', 'Products', ...
    'Nodes', 'SolverStatus', ...
    'ObjectiveErrorBound', 'IntegerObjective'});
metadata = struct('matlab', version, 'host', getenv('HOSTNAME'), ...
    'job', getenv('SLURM_JOB_ID'), 'array', getenv('SLURM_ARRAY_JOB_ID'), ...
    'total_limit', total_limit, 'p1_limit', p1_limit, ...
    'input_file', cases.File(k), 'input_sha256', cases.SHA256(k), ...
    'construction', 'continuous-bound', 'phase2', 'direct-miqp', ...
    'gurobi_interface', which('gurobi'));
if isfield(answer, 'version')
    metadata.gurobi = answer.version;
end
save([file, '.mat'], 'row', 'r', 'd', 'model', 'answer', 'checks', 'pairs', ...
    'failure', 'params', 'options', 'metadata', '-v7.3');
writetable(row, [file, '.csv']);
disp(row);
assert(status ~= "ERROR", 'Method failed; evidence saved.');
end
