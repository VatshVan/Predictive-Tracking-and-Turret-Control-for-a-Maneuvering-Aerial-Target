%PARITY_CHECK Replay logged Python controller inputs through turret_controller.m.
%   Run python/scripts/export_parity.py first. For every parity/*.csv the recorded
%   inputs are fed to turret_controller.m and its outputs are compared with the
%   Python outputs (open-loop replay, so any difference shows up directly).

here = fileparts(mfilename('fullpath'));
files = dir(fullfile(here, 'parity', '*.csv'));
if isempty(files)
    error('No parity files: run python/scripts/export_parity.py first.');
end
all_ok = true;
fprintf('%-16s %8s %14s %14s %10s\n', 'case', 'calls', 'max|dV| [V]', 'max|dV|/|V|', 'fire diff');
for f = files'
    M = readmatrix(fullfile(f.folder, f.name));
    clear turret_controller
    n = size(M, 1);
    V = zeros(n, 1); fire = zeros(n, 1);
    for k = 1:n
        [V(k), fire(k)] = turret_controller(M(k,1), M(k,2), M(k,3), M(k,4), M(k,5), M(k,6));
    end
    dV = abs(V - M(:,7));
    rel = max(dV ./ max(abs(M(:,7)), 1));
    nfire = sum(fire ~= M(:,8));
    ok = max(dV) < 1e-6 && nfire == 0;
    all_ok = all_ok && ok;
    fprintf('%-16s %8d %14.3e %14.3e %10d %s\n', f.name, n, max(dV), rel, nfire, string(ok));
end
if all_ok
    fprintf('PARITY OK: MATLAB and Python controllers agree.\n');
else
    fprintf('PARITY FAILED\n');
end
