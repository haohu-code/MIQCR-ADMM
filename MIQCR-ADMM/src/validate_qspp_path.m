function validate_qspp_path(s, x)
% Check binary flow feasibility and a single simple source-to-sink path.
assert(all(ismember(x, [0, 1])) && norm(s.incidence*x-s.rhs, Inf) < 1e-8);
used = find(x);
pos = s.source;
seen = pos;

% Follow the selected arcs, rejecting branches and repeated vertices.
while pos ~= s.sink
    arcs = used(s.incidence(pos, used) == 1);
    assert(numel(arcs) == 1, 'Not a simple path');
    e = arcs(1);
    used(used == e) = [];
    pos = find(s.incidence(:, e) == -1);
    assert(~ismember(pos, seen));
    seen(end+1) = pos;
end
assert(isempty(used), 'Extra arcs outside path');
end
