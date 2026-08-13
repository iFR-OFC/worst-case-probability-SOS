function plot_ex2_sos_ambiguity(create_tikz)

%% PARSE INPUT ARGUMENTS
tikz = nargin > 0 && create_tikz;

%% LOAD DATA
data = load("data/ex2_sos_ambiguity_results.mat");
results = data.results;

%% PLOTTING PARAMETERS
% color scheme (greyscaled for publication)
color1 = '#242424'; % blue

markers = {'o', 's', 'd', '^', 'v'};
line_styles = {'-', '--', ':', '-.', '-'};

%% CREATE FIGURE
fighandle = figure();

%% PLOT UPPER BOUND 
for i=1:length(results.list_ambiguity)
    plot(results.list_degree, results.upper_bound(:,i), line_styles{i}, 'LineWidth', 1, ...
        'Marker', markers{i}, 'MarkerSize', 8, 'Color', color1,...
        'DisplayName', sprintf('\\alpha = %.3f', results.list_ambiguity(i)));
    hold on
end
ylabel('P_d*', 'FontSize', 12, 'Color', 'k');
ylim([0.5, 1]);
grid on;

% global setting
xlabel('d', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 11);

%% EXPORT TO TIKZ
if tikz == 1
    matlab2tikz('figures/ex2_stats_sos_ambiguity.tex', 'figurehandle', fighandle, 'standalone', true);
end

end