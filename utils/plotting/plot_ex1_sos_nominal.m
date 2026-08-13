function plot_ex1_sos_nominal(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex1_sos_nominal_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; 
color2 = '#4f4f4f'; 
color3 = '#797979'; 
color4 = '#ababab';

%% CREATE FIGURE
fighandle = figure();

% get data from struct
list_degree = results.list_degree;
stats_list  = results.status;
value_list  = results.upper_bound;
built_time  = results.built_time; 
solve_time  = results.solve_time;
prob_list   = results.lower_bound;

%% Plot upper bound and potential reference (left y-axis)
yyaxis left;
plot(list_degree, value_list, '-', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color1);
hold on
h = scatter(list_degree(stats_list), value_list(stats_list), 80, 'o', 'MarkerEdgeColor', color1);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~stats_list), value_list(~stats_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('P_d*', 'FontSize', 12, 'Color', 'k');
ylim([0, 1]);
grid on;

% plot the potential upper bound through a test case
yyaxis left
plot(list_degree, max(prob_list)*ones(length(list_degree),1), '-', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color4);

%% Plot computation times (right y-axis)
yyaxis right;
semilogy(list_degree, built_time, '-.', 'LineWidth', 1, 'Color', color3);
hold on
semilogy(list_degree, solve_time, '--', 'LineWidth', 1, 'Color', color3);

h = scatter(list_degree(stats_list), built_time(stats_list), 80, 'square', 'MarkerEdgeColor', color3);

h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~stats_list), built_time(~stats_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

h = scatter(list_degree(stats_list), solve_time(stats_list), 80,'square', 'MarkerEdgeColor', color3);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h = scatter(list_degree(~stats_list), solve_time(~stats_list), 80, 'diamond', 'filled', ...
    'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('time (s)', 'FontSize', 12, 'Color', 'k');

% Dummy object for the legend
if any(~stats_list)
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
    matlab2tikz('figures/ex1_stats_sos_nominal.tex', 'figurehandle', fighandle, 'standalone', true);
end

end