function run_all_experiments()
% Execute all experiments and save results

addpath("experiments/")

% add required paths
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

% create data directory if it doesn't exist
data_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

% get experiments to run
root_dir = fileparts(fileparts(mfilename('fullpath')));
exp_dir = fullfile(root_dir, 'experiments');
% check if directory exists
if ~exist(exp_dir, 'dir')
    error('Experiments directory not found: %s', exp_dir);
end

% get all files
files = dir(exp_dir);
files = files(~[files.isdir]);
experiments = {files.name};
experiments = strrep(experiments, '.m', '');
num_experiments = length(experiments);

% display start message
fprintf('========================================\n');
fprintf('Starting all experiments...\n');
fprintf('========================================\n\n');

success_count = 0;

% run each experiment
for i = 1:num_experiments
    % remove the .m termination
    exp_name = experiments{i};

    fprintf('Running experiment %d/%d:\n', i, num_experiments);
    fprintf('  Function: %s\n', exp_name);
    
    try
        % run the experiment
        tic;
        feval(exp_name);
        elapsed_time = toc;
        
        fprintf('  ✓ Success! (%.2f seconds)\n', elapsed_time);
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
    clearvars -except exp_name experiments i num_experiments success_count
end

% display summary
fprintf('========================================\n');
fprintf('Experiment Summary\n');
fprintf('========================================\n');
fprintf('Total: %d\n', num_experiments);
fprintf('Successful: %d\n', success_count);
fprintf('Failed: %d\n', num_experiments - success_count);
fprintf('========================================\n');

end