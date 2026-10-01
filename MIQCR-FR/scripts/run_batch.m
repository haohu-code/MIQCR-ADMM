function run_batch(tasks, folder, total_limit, p1_limit)
% Run selected tasks sequentially; preserve failures and continue with the rest.
% Example: run_batch(1:18, 'outputs/example', 300, 120).
root = setup;
if nargin < 1 || isempty(tasks)
    cases = readtable(fullfile(root, 'data/manifest.csv'), ...
        'TextType', 'string', 'Delimiter', ',', 'NumHeaderLines', 0);
    tasks = 1:3*height(cases);
end
if nargin < 2
    folder = 'outputs/batch';
end
if nargin < 3
    total_limit = 3600;
end
if nargin < 4
    p1_limit = 900;
end
for task = tasks(:)'
    try
        run_experiment(task, total_limit, p1_limit, folder);
    catch err
        fprintf(2, 'Task %d failed or already exists: %s\n', task, err.message);
    end
end
% The collector identifies missing, failed, and inconsistent records.
end
