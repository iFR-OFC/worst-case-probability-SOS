function plot_ex2_sos_multiple(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex2_sos_multiple_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; 

line_styles = {'-', '--', ':', '-.', '-'};
markers = {'o', 's', 'd', '^', 'v'};

%% CREATE FIGURE
fighandle = figure();

%% PLOT UPPER BOUND
for i=1:length(results.list_Kcenter)
    plot(results.list_degree, results.upper_bound(:,i), line_styles{i}, 'LineWidth', 1, ...
        'MarkerSize', 8, 'Marker', markers{i}, 'Color', color1, ...
        'DisplayName', sprintf('K_%d', i));
    hold on
end
ylabel('Upper bound', 'FontSize', 12, 'Color', 'k');
ylim([0.5, 1]);
grid on;

% global setting
xlabel('Relaxation order d', 'FontSize', 12);
%legend('Upper bound', 'Build time', 'Solve time', ' Solver got issues');
legend('Orientation', 'vertical');

%% EXPORT TO TIKZ
if tikz == 1
    matlab2tikz('figures/ex2_stats_sos_multiple.tex', 'figurehandle', fighandle, 'standalone', true);
end

end