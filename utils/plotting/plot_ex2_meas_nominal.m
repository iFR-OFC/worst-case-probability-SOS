function plot_ex2_meas_nominal(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex2_meas_nominal_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; 
color2 = '#4f4f4f'; 

%% CREATE FIGURE
fighandle = figure();

%% PLOT UPPER BOUND (left y-axis)
yyaxis left;
plot(results.list_degree, results.upper_bound, '-o', 'Color', color1, 'LineWidth', 1, 'MarkerSize', 8);
ylabel('P_d*', 'FontSize', 12, 'Color', color1);
ylim([0, 1]);
grid on;

%% PLOT COMPUTATION TIMES (right y-axis)
yyaxis right;
semilogy(results.list_degree, results.built_time, '-.s', 'Color', color2, 'LineWidth', 1, 'MarkerSize', 8);
hold on
semilogy(results.list_degree, results.solve_time, '--s', 'Color', color2, 'LineWidth', 1, 'MarkerSize', 8);
ylabel('time (s)', 'FontSize', 12, 'Color', color2);

% global setting
xlabel('d', 'FontSize', 12);
legend('Upper bound', 'Build time', 'Solve time');

%% EXPORT TO TIKZ
if tikz == 1
    matlab2tikz('figures/ex2_stats_meas_nominal.tex', 'figurehandle', fighandle, 'standalone', true);
end

end