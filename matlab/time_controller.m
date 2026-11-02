%TIME_CONTROLLER Measure turret_controller.m execution time per call (scenario B, 1 kHz).
%   Sets the upper bound on the controller call rate quoted in the README.

f_ctrl = 1000;
o = run_engagement('B', f_ctrl);          % closed-loop run to obtain realistic inputs
clear turret_controller
n = numel(o.t);
dt_call = zeros(n, 1);
for k = 1:n
    [r_m, th_m, valid] = target_motion_generator(o.t(k), 'B');
    tic;
    turret_controller(o.t(k), r_m, th_m, valid, o.th(k), o.thd(k));
    dt_call(k) = toc;
end
fprintf('turret_controller.m: mean %.3f ms, median %.3f ms, p99 %.3f ms, max %.3f ms per call\n', ...
        1e3*mean(dt_call), 1e3*median(dt_call), 1e3*prctile(dt_call, 99), 1e3*max(dt_call));
