%% Basic capacity consensus
%
% This is the basic capacity consensus that is used as a baseline; it
% does not support either delays nor has a finite termination condition -
% meaning: that is guaranteed to terminate as our main algorithm.
%
% Authors:
%
%  - Andreas A. Grammenos (ag926@cl.cam.ac.uk)
%  - Themistoklis Charalambous (themistoklis.charalambous@aalto.fi)
%
% License: GPLv3
%
%

%% Initialisation
%
clear; clc; close all;

% for reproducibility, fix the seed
rng("default")

%% Configuration
%

% nodes to test for
nodes_to_test = [20, 50, 100, 200, 300, 600];
node_test_len = length(nodes_to_test);
node_len_array = 1:node_test_len;

% Number of iterations
max_iter = 30;

% plotted for trial
plotted_for_trial = 0;

% number of trials
trials = 30;

% enable extensive plots
extensive_plots = 1;

% multiplier from sec to ns & μs
ns_mul = 1e+9;
micro_mul = 1e+6;

% the epsilon to check if we converged
e = 1e-5;

% the minimum workload value
min_workload = 1;
% the maximum workload value
max_workload = 2;

% timers
total_time_global = zeros(trials, node_test_len);
iter_time_global = zeros(trials, node_test_len, max_iter);

% converge statistics
cov_min_global = zeros(node_test_len, trials);
cov_max_global = zeros(node_test_len, trials);
cov_mean_global = zeros(node_test_len, trials);
cov_win_global = zeros(node_test_len, trials);

% setup variables
params.type = "basic";
params.pflag = 0;
params = setup_vars(params);

%% Run the trials
%

for t=1:trials
  fprintf("\n** Running trial %d\n", t);
  % run for the specified configuration
  for i=node_len_array

    % number of nodes to consider
    nodes = nodes_to_test(i);
    
    % note the converged time for each node
    node_converged_at = zeros(nodes, 1);

    fprintf("\t-- INFO: Testing for %d nodes in the graph\n", nodes);

    % generate the graph
    [~, ~, n, P] = gen_graph(nodes); 

    % generate the workload to be shared among users
    x0 = gen_workload(min_workload, max_workload, nodes);

    % create the node capacities

    y1= 3 * ones(ceil(nodes / 2), 1);
    y2= 5 * ones(floor(nodes / 2), 1);

    y = [y1; y2];

    % generate the "booked" utilization
    u = gen_utilisation(nodes);

    x = x0 + u;
    f = [y1; y2];

    archive_z = [0 x'; zeros(max_iter, nodes+1)];

    % create the time for the iterations to be saved
    iter_run_time = zeros(max_iter, 1);

    % -- Run the consensus algorithm -- %

    % Perform the basic consensus algorithm (standard)
    total_time = tic;
    for k=1:max_iter 
      % start the tag
      iter_time = tic;
      x = P*x;
      y = P*y;
      z = (x./y).*f-u;
      % finish the tag
      iter_run_time(k) = micro_mul * toc(iter_time);
      
      % record the numbers
      archive_z(k+1, :) = [k z'];
      
      % check about convergence - leave first row as that's 
      % the iteration number.
      delta = abs(archive_z(k+1, 2:end) - archive_z(k, 2:end))';
      % check which converged
      converged = find(delta < e);
      % check which have not converged yet
      non_converged_yet = find(node_converged_at == 0);
      
      % only check for a match if we have an index that has 
      % been converged but also we have non-converted indices as well
      if ~isempty(converged) && ~isempty(non_converged_yet)
        % this part is quite tricky - what this does is takes the 
        % intersection of the converted with the non_converted_yet indices
        % and finds which ones are present in both!
        %
        % However, we are not interested in the actual indices, rather than
        % the values that they hold, which are stored in the idcs variable.
        [idcs, ~] = intersect(non_converged_yet, converged, 'stable');
        % then using logical indexing we update the converge time for these
        % nodes to be equal to the current iteration.
        node_converged_at(idcs) = k;
      end
      
    end
    total_time = micro_mul * toc(total_time);

    % convergence statistics
    cov_min = min(node_converged_at);
    cov_max = max(node_converged_at);
    cov_mean = mean(node_converged_at);
    cov_win = cov_max - cov_min;
    fprintf("\t== RUN INFO: mean conv. iter: %.2f, min: %d, max: %d, window: %d\n", ...
      cov_mean, cov_min, cov_max, cov_win);
    
    % store the aggregate statistics
    cov_min_global(i, t) = cov_min;
    cov_max_global(i, t) = cov_max;
    cov_mean_global(i, t) = cov_mean;
    cov_win_global(i, t) = cov_win;
    
    % -- Plot the results -- %

    if (extensive_plots == 1 || nodes <= 10) && plotted_for_trial == 0
      fig = figure;
      hold on;
        plot(0:max_iter, archive_z(:, 2:nodes+1));
      hold off
      title("Values at each node vs Number of iterations", ...
         "Interpreter", "Latex", "FontName", "Times New Roman");
      xlabel("Number of iterations", ...
        "Interpreter", "Latex", "FontName", "Times New Roman");
      ylabel("Value at each node", ...
        "Interpreter", "Latex", "FontName", "Times New Roman");
      % print the figure 
      st = sprintf("n_%d_iter_%d_values_for_each_node", nodes, k);
      print_fig(fig, st, params);
      
      % plot the execution time
      fig = figure;
      plot(1:max_iter, iter_run_time);
      title(sprintf("Time per iteration with %d nodes", nodes), ...
         "Interpreter", "Latex", "FontName", "Times New Roman");
      xlabel("Number of iterations", ...
         "Interpreter", "Latex", "FontName", "Times New Roman");
      ylabel("Time ($\mu s$)", ...
         "Interpreter", "Latex", "FontName", "Times New Roman");
      % print the figure 
      st = sprintf("n_%d_iter_%d_iter_execution_time", nodes, k);
      print_fig(fig, st, params);
    end

    % save global stats
    iter_time_global(t, i, :) = iter_run_time;
    total_time_global(t, i) = total_time;

    fprintf("\t-- INFO: Finished %d nodes with exec %.3f μs and mean iter: %.3f μs\n", ...
      nodes, total_time, mean(iter_run_time));
  end
  % raise it to avoid printing for the rest of the trials
  plotted_for_trial = 1;
  fprintf("** Finished trial %d\n", t);
end
fprintf("\n");

%% Final (aggregate) plots
%

% plot the mean execution time per node batch for all trials
fig = figure;
std_ttg = std(total_time_global) / sqrt(trials);
avg_ttg = mean(total_time_global);
errorbar(node_len_array, avg_ttg, std_ttg, '-*', 'LineWidth', 2);
title(sprintf("Average total time spend with %d iterations over %d trials", ... 
  max_iter, trials), "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(node_len_array)
xticklabels(num2cell(nodes_to_test))
xlabel("Nodes", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
ylabel("Time ($\mu s$)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
% print the figure 
st = sprintf("n_%d_iter_%d_trials_%d_total_exec_time", nodes, k, trials);
print_fig(fig, st, params);

% plot the mean iteration execution time for node batch for all trials
fig = figure;
hold on;
wlegs = cell(1, node_test_len);
for i=node_len_array
  cur_iter = squeeze(iter_time_global(i, :, :));
  std_cur = std(cur_iter) / sqrt(trials);
  mean_cur = mean(cur_iter);
  % plot looks cleaner for this, but leaving errorbar as well.
  % errorbar(mean_cur, std_cur, 'LineWidth', 2);
  plot(mean_cur, 'LineWidth', 2);
  wlegs{i} = sprintf("N=%d", nodes_to_test(i));
end
title("Average iteration time", ...
  "Interpreter", "Latex", "FontName", "Times New Roman")
hold off;
xlabel("Iterations", ...
   "Interpreter", "Latex", "FontName", "Times New Roman");
ylabel("Time ($\mu s$)", ...
   "Interpreter", "Latex", "FontName", "Times New Roman");
legend(wlegs);
% print the figure 
st = sprintf("n_%d_iter_%d_trials_%d_avg_iter_exec_time", nodes, k, trials);
print_fig(fig, st, params);

% plot the mean convergence for each node batch
fig = figure;
hold on;
  % -- MIN CONVERGENCE -- %
  min_avg = mean(cov_min_global, 2);
  min_std = std(cov_min_global, 0, 2) / sqrt(trials);
  errorbar(node_len_array, min_avg, min_std, '-*', 'LineWidth', 2);
  % -- MAX CONVERGENCE -- %
  max_avg = mean(cov_max_global, 2);
  max_std = std(cov_max_global, 0, 2) / sqrt(trials);
  errorbar(node_len_array, max_avg, max_std, '-*', 'LineWidth', 2);
  % -- MEAN CONVERGENCE -- %
  avg_avg = mean(cov_mean_global, 2);
  avg_std = std(cov_mean_global, 0, 2) / sqrt(trials);
  errorbar(node_len_array, avg_avg, avg_std, '-*', 'LineWidth', 2);
  % -- WINDOW OF CONVERGENCE -- %
  win_avg = mean(cov_win_global, 2);
  win_std = std(cov_win_global, 0, 2) / sqrt(trials);
  errorbar(node_len_array, win_avg, win_std, '-*', 'LineWidth', 2);
hold off;

xticks(node_len_array)
xticklabels(num2cell(nodes_to_test))
xlabel("Nodes", ...
   "Interpreter", "Latex", "FontName", "Times New Roman");
ylabel("Iteration", ...
   "Interpreter", "Latex", "FontName", "Times New Roman");
title("Convergence statistics", ...
   "Interpreter", "Latex", "FontName", "Times New Roman");
legend("min", "max", "mean", "window");
% print the figure 
st = sprintf("n_%d_iter_%d_trials_%d_agg_cov_statistics", nodes, k, trials);
print_fig(fig, st, params);

