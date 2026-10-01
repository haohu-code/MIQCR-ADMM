function root = setup()
% Initialize the release independently of the enclosing research checkout.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'), fullfile(root, 'scripts'));

% Use the configured Gurobi installation, with the existing macOS fallback.
home = getenv('GUROBI_HOME');
if ~isempty(home)
    addpath(fullfile(home, 'matlab'));
end
if exist('gurobi', 'file') == 0 && ismac && ...
        isfolder('/Library/gurobi1100/macos_universal2/matlab')
    addpath('/Library/gurobi1100/macos_universal2/matlab');
end
assert(exist('gurobi', 'file') ~= 0, ...
    'Set GUROBI_HOME or add the Gurobi MATLAB interface to the path.');
end
