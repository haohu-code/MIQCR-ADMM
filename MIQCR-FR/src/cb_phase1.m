function [d, seconds] = cb_phase1(p, folder, limit, hard_limit)
% Run the external MIQCR-CB Phase 1 and read its existing JSON export.

clock = tic; % Include input preparation and export validation in P1.

% Locate the dependency and write the original binary quadratic problem.
root = setup();
external = getenv('MIQCR_EXTERNAL');
if isempty(external)
    external = fullfile(fileparts(root), 'external');
end
dependency = fullfile(external, 'MIQCR-CB');
if ~isfolder(dependency)
    dependency = fullfile(external, 'smiqp');
end
alg = fullfile(dependency, 'source', 'Smiqp-1.0', 'src', 'Alg');
if ~isfolder(folder)
    mkdir(folder);
end
file = fullfile(folder, [p.name, '.dat']);
export = fullfile(folder, [p.name, '.json']);
checkpoint = [export, '.checkpoint.json'];
log = fullfile(folder, [p.name, '.log']);
assert(~isfile(file) && ~isfile(export) && ~isfile(log) && ~isfile(checkpoint), ...
    'Refuse to overwrite CB artifacts; use a fresh output directory.');
write_smiqp(p, file);

% Enforce P1 externally; only complete atomic disk checkpoints survive a kill.
executable = './smiqp';
limit_env = ['SMIQP_CHECKPOINT_PATH=',quote(checkpoint),' '];
if nargin >= 3
    validateattributes(limit, {'numeric'}, {'scalar','finite','positive'});
    if nargin >= 4
        limit = min(limit, hard_limit);
    end
    remaining = limit-toc(clock);
    if remaining <= 0
        d = [];
        seconds = toc(clock);
        return;
    end
    limit_env = [limit_env,sprintf('SMIQP_CB_WALL_LIMIT=%.17g ',remaining)];
    executable = sprintf('python3 %s --hard %g ./smiqp', ...
        quote(fullfile(root,'scripts','run_with_limit.py')),remaining);
end
command = ['cd ', quote(alg), ' && env SMIQP_PHASE1_ONLY=1 ', ...
    'SMIQP_PHASE2_BETA_TOLERANCE=0.0001 SMIQP_EXPORT_PHASE2=', quote(export), ...
    ' ', limit_env, executable, ' ', quote(file), ' > ', quote(log), ' 2>&1'];
status = system(command);
seconds = toc(clock);
% Only an atomically completed export/checkpoint is accepted after a hard stop.
timed_out = nargin >= 3 && any(status == [124, 137]);
assert(status == 0 || timed_out, 'MIQCR-CB failed; inspect its log.');
logtext = fileread(log);
completed = contains(logtext, 'Phase-1-only run complete.');
if nargin >= 3
    assert(contains(logtext, 'CB graceful wall limit:'), ...
        'Rebuild MIQCR-CB: executable lacks the graceful stopping interface.');
end
soft_stop = contains(logtext, 'CB graceful time limit reached; center evaluated.');
used_checkpoint = false;
if timed_out && (~isfile(export) || ~completed)
    if ~isfile(checkpoint)
        d = [];
        return;
    end
    d = read_cb_export(p,checkpoint);
    used_checkpoint = true;
else
    assert(isfile(export) && completed, 'MIQCR-CB did not complete its export.');
    d = read_cb_export(p,export);
end
d.used_disk_checkpoint = used_checkpoint;
d.process_status = status;
d.oracle_fallback = contains(logtext, 'CB oracle fallback: restored last completed evaluation');
d.reduced_accuracy_oracles = numel(strfind(logtext, ...
    'CB SDP oracle: reduced accuracy (status 3).'));
d.timed_out = timed_out || soft_stop;
d.soft_time_limit = soft_stop;
d.bundle_termination = regexp(logtext, ...
    'termination status: ([^\r\n]+)', 'tokens', 'once');
if ~isempty(d.bundle_termination)
    d.bundle_termination = d.bundle_termination{1};
else
    d.bundle_termination = 'see log';
end
if used_checkpoint
    d.termination = 'TIME_LIMIT_WITH_CHECKPOINT';
elseif d.oracle_fallback
    d.termination = 'ORACLE_FAILURE_WITH_SNAPSHOT';
elseif soft_stop
    d.termination = 'SOFT_TIME_LIMIT_WITH_EXPORT';
elseif timed_out
    d.termination = 'TIME_LIMIT_WITH_EXPORT';
else
    d.termination = d.bundle_termination;
end
seconds = toc(clock); % Include export validation in the Phase-1 budget.
end

function s = quote(s)
% Preserve literal file paths in the shell command.
a = char(39);
s = [a, strrep(s, a, char([39, 34, 39, 34, 39])), a];
end
