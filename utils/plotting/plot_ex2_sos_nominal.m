function plot_ex2_sos_nominal(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex2_sos_nominal_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; 
color2 = '#4f4f4f'; 
color3 = '#797979'; 

%% CREATE FIGURE
fighandle = figure();

list_degree     = results.list_degree;
value_list      = results.upper_bound;
status_list     = results.status;
built_time_list = results.built_time;
solve_time_list = results.solve_time;

%% PLOT UPPER BOUND (left y-axis)
yyaxis left;
plot(list_degree, value_list, '-', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color1);
hold on
h = scatter(list_degree(status_list), value_list(status_list), 80, 'o', 'MarkerEdgeColor', color1);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~status_list), value_list(~status_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('P_d*', 'FontSize', 12, 'Color', 'k');
ylim([0, 1]);
grid on;

%% PLOT COMPUTATION TIMES (right y-axis)
yyaxis right;
semilogy(list_degree, built_time_list, '-.', 'LineWidth', 1, 'Color', color3);
hold on
semilogy(list_degree, solve_time_list, '--', 'LineWidth', 1, 'Color', color3);

h = scatter(list_degree(status_list), built_time_list(status_list), 80, 'square', 'MarkerEdgeColor', color3);

h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~status_list), built_time_list(~status_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

h = scatter(list_degree(status_list), solve_time_list(status_list), 80,'square', 'MarkerEdgeColor', color3);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~status_list), solve_time_list(~status_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('Time (s)', 'FontSize', 12, 'Color', 'k');

% Dummy object for the legend
if any(~status_list)
    h = plot([NaN NaN], [NaN NaN], ...
        'LineStyle', 'none', ...
        'Marker', 'd', ...
        'MarkerFaceColor', color2, ...
        'MarkerEdgeColor', color2, ...
        'MarkerSize', 8);
    % Add a second marker (blue square) to the same object
    h.MarkerIndices = [1 2];   % MATLAB doesn't allow different marker types/colors

    legend('Upper bound', 'Build time', 'Solve time', ' Solver got issues');
else
    legend('Upper bound', 'Build time', 'Solve time');
end

ax = gca;
ax.YAxis(1).Color = color1;
ax.YAxis(2).Color = color3;
ax.YAxis(2).Scale = 'log';

% global setting
xlabel('Relaxation order d', 'FontSize', 12);
legend('Orientation', 'horizontal', 'Location', 'northoutside');

%% EXPORT TO TIKZ
if tikz == 1
    matlab2tikz('figures/ex2_stats_sos_nominal.tex', 'figurehandle', fighandle, 'standalone', true);
end

end