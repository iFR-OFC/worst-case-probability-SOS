function generate_all_figures()
% Generate all figures from saved results

addpath("utils/plotting/")

% add required paths
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));

% define directories
fig_dir = fullfile(root_dir, 'figures');

% create figures directory if it doesn't exist
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

plotgen_dir = fullfile(root_dir, 'utils/plotting');

% get all files
files = dir(plotgen_dir);
files = files(~[files.isdir]);
plotgen = {files.name};
plotgen = strrep(plotgen, '.m', '');
num_plotgen = length(plotgen);

% display start message
fprintf('========================================\n');
fprintf('Generating all figures...\n');
fprintf('========================================\n\n');

success_count = 0;

for i = 1:num_plotgen
    plot_func = plotgen{i};
    
    fprintf('Generating figure %d/%d: \n', i, num_plotgen);
    fprintf('  Plot function: %s\n', plot_func);
    
    try
        % call the plotting function
        tic
        feval(plot_func);
        elapsed_time = toc;
        
        fprintf('  ✓ Success! (%.2f seconds)\n', elapsed_time);
        fprintf('  Figures .tex saved to: %s/\n\n', fig_dir);
        success_count = success_count + 1;
        
    catch ME
        fprintf('  ✗ Failed!\n');
        fprintf('  Error: %s\n', ME.message);
        fprintf('  Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('\n');
    end
    clearvars -except num_plotgen success_count i fig_dir plotgen
end

% display final summary
fprintf('========================================\n');
fprintf('Figure Generation Summary\n');
fprintf('========================================\n');
fprintf('Total tasks: %d\n', num_plotgen);
fprintf('Successful: %d\n', success_count);
fprintf('Failed: %d\n', num_plotgen - success_count);
fprintf('========================================\n');

end
