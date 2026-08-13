function plot_ex1_meas_nominal(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex1_meas_nominal_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; 
color2 = '#4f4f4f'; 
color3 = '#797979'; 
color4 = '#ababab';

%% CREATE FIGURE
fighandle = figure();

list_degree = results.list_degree;
has_failures = any(~results.status);

%% Plot upper bound and potential reference (left y-axis)
yyaxis left;
hold on
% Upper bound curve
plot(list_degree, results.upper_bound, '-', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color1);
% Scatter points for upper bound
h = scatter(list_degree(results.status), results.upper_bound(results.status), 80, 'o', 'MarkerEdgeColor', color1);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~results.status), results.upper_bound(~results.status), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% labels
ylabel('P_d*', 'FontSize', 12, 'Color', 'k');
ylim([0, 1]);
grid on;

% Reference line: lower bound obtained via Monte Carlo simulations
yyaxis left
plot(list_degree, max(results.lower_bound)*ones(length(list_degree),1), '-', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color4);

%% Plot computation times (right y-axis)
yyaxis right;
hold on
semilogy(list_degree, results.built_time, '-.', 'LineWidth', 1, 'Color', color3);
semilogy(list_degree, results.solve_time, '--', 'LineWidth', 1, 'Color', color3);

% Scatter points for build time
h = scatter(list_degree(results.status), results.built_time(results.status), 80, 'square', 'MarkerEdgeColor', color3);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~results.status), results.built_time(~results.status), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Scatter points for solve time
h = scatter(list_degree(results.status), results.solve_time(results.status), 80,'square', 'MarkerEdgeColor', color3);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~results.status), results.solve_time(~results.status), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('time (s)', 'FontSize', 12, 'Color', 'k');

% Dummy object for the legend
if has_failures
    h = plot([NaN NaN], [NaN NaN], ...
        'LineStyle', 'none', ...
        'Marker', 'd', ...
        'MarkerFaceColor', color2, ...
        'MarkerEdgeColor', color2, ...
        'MarkerSize', 8);
    % Add a second marker (blue square) to the same object
    h.MarkerIndices = [1 2];   % MATLAB doesn't allow different marker types/colors

    legend('Upper bound', 'Arbitrary initial measure', 'Build time', 'Solve time' ,' Solver got issues');
else
    legend('Upper bound', 'Arbitrary initial measure', 'Build time', 'Solve time');
end

%% Configure axes
ax = gca;
ax.YAxis(1).Color = color1;
ax.YAxis(2).Color = color3;
ax.YAxis(2).Scale = 'log';

% global setting
xlabel('d', 'FontSize', 12);

%% Export to TikZ if requested
if tikz == 1
    matlab2tikz('figures/ex1_stats_meas_nominal.tex', 'figurehandle', fighandle, 'standalone', true);
end

end