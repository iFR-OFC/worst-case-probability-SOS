function plot_ex1_sos_ambiguity(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex1_sos_ambiguity_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1    = '#242424';  
markers   = {'o', 's', 'd', '^', 'v'};
linestyle = {'-', '--', ':', '-.', '-'};

%% CREATE FIGURE
fighandle = figure();

% Plot computation time
for i=1:length(results.list_ambiguity)
    plot(results.list_degree, results.upper_bound(:,i), ...
        linestyle{i}, ...
        'LineWidth', 1, ...
        'Marker', markers{i}, ...
        'MarkerSize', 8, ...
        'Color', color1,...
        'DisplayName', sprintf('\\alpha = %.3f', results.list_ambiguity(i)));
    hold on
end
% label
ylabel('P_d*', 'FontSize', 12, 'Color', 'k');
xlabel('d', 'FontSize', 12);

% limits
ylim([0, 1]);
xlim([1,8])

% global setting
legend('Location', 'best', 'FontSize', 11);
grid on;

%% Export to TikZ if requested
if tikz == 1
    matlab2tikz('figures/ex1_stats_sos_ambiguity.tex', 'figurehandle', fighandle, 'standalone', true);
end

end