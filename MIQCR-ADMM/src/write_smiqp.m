function write_smiqp(p, file)
% Write the binary quadratic problem in SMIQP's zero-based input format.
n = numel(p.c);
fid = fopen(file, 'w');
assert(fid >= 0);

% Dimensions and binary upper/lower bounds.
fprintf(fid, '%d\n%d\n%d\n%d\n0\n0\nu\n', n, n, size(p.G, 1), size(p.A, 1));
fprintf(fid, '%.0f ', ones(1, n));
fprintf(fid, '\nl\n');
fprintf(fid, '%.0f ', zeros(1, n));
fprintf(fid, '\n');

% Objective and original linear equalities/inequalities.
matrix(fid, 'q', p.Q);
vector(fid, 'c', p.c);
if ~isempty(p.h)
    matrix(fid, 'a', p.G);
    vector(fid, 'b', p.h);
end
if ~isempty(p.b)
    matrix(fid, 'd', p.A);
    vector(fid, 'e', p.b);
end
fclose(fid);
end

function matrix(fid, name, A)
[i, j, v] = find(A);
fprintf(fid, '%s\n%d\n', name, numel(v));
if ~isempty(v)
    fprintf(fid, '%d %d %.17g\n', [i(:)-1, j(:)-1, v(:)]');
end
end

function vector(fid, name, v)
v = full(v);
i = find(v);
fprintf(fid, '%s\n%d\n', name, numel(i));
if ~isempty(i)
    fprintf(fid, '%d %.17g\n', [i(:)-1, v(i(:))]');
end
end
