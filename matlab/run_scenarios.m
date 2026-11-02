%RUN_SCENARIOS Simulate the standardized and stress scenarios and plot each engagement.
%   Edit f_ctrl to change the controller call rate (20 Hz .. see README).

f_ctrl = 100;
scenarios = {'A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6'};
here = fileparts(mfilename('fullpath'));
fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('%-4s %8s %12s %12s %10s %8s\n', 'case', 'hits', 'first hit', 'point RMS', 'in R_eng', 'V peak');
for c = 1:numel(scenarios)
    sc = scenarios{c};
    o = run_engagement(sc, f_ctrl);
    m = o.metrics;
    fprintf('%-4s %4d/%-3d %10.3f s %10.2f m %10.3f %6.0f V\n', sc, m.n_hits, m.n_pulses, ...
            m.t_first_hit, m.point_rms_m, m.frac_in_Reng, m.V_peak);

    fig = figure('Visible', 'off', 'Position', [100 100 560 760]);
    tiledlayout(4, 1, 'TileSpacing', 'compact');
    nexttile; hold on; grid on
    plot(o.t, rad2deg(o.th_true), 'k', 'LineWidth', 1.0);
    plot(o.t, rad2deg(o.th_hat), '--', 'Color', [0.16 0.47 0.84]);
    plot(o.t, rad2deg(o.r_ref), ':', 'Color', [0.92 0.41 0.20], 'LineWidth', 1.2);
    plot(o.t, rad2deg(o.th), 'Color', [0.11 0.69 0.48], 'LineWidth', 1.4);
    ylabel('elevation [deg]');
    legend({'true', 'estimate', 'reference \theta(t+t_p)', 'turret'}, 'Location', 'southeast');
    title(sprintf('Scenario %s: first hit %.2f s, %d/%d pulses hit', sc, m.t_first_hit, m.n_hits, m.n_pulses));
    nexttile; hold on; grid on
    err = o.r_true .* sin(o.th - o.th_true);
    plot(o.t, err, 'Color', [0.11 0.69 0.48]);
    yline([-2 2], '--k');
    if ~isempty(o.pulses)
        k = round(o.pulses(:,1) * f_ctrl) + 1;
        plot(o.pulses(:,1), err(k), 'ko', 'MarkerFaceColor', 'k');
    end
    ylim([-6 6]); ylabel('cross-range [m]');
    nexttile; grid on
    plot(o.t, o.V, 'Color', [0.16 0.47 0.84]); ylabel('V_m [V]'); ylim([-650 650]);
    nexttile; grid on
    plot(o.t, o.p_hit, 'Color', [0.92 0.41 0.20]); ylabel('P_{hit}'); ylim([-0.05 1.05]);
    xlabel('time [s]');
    exportgraphics(fig, fullfile(fig_dir, sprintf('scenario_%s.png', sc)), 'Resolution', 200);
    close(fig);
end
