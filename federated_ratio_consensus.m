%% Federated capacity ratio consensus with finite time termination
%
% This is the main stub that performs the experiments in our paper. We test
% across a variety of parameters; concretely we define a parameter tuple as
% the node count, delay and run it across a number of trials - in the paper
% we used 10 trials.
%
% We show, that every node in the network converges within multiples of 
% (1+τ)D, where τ is the delay and D is the diameter (i.e.: distance) of 
% the network. Our experiments can be considered large scale as they scale 
% from 20 to 600 nodes in a single network and each node can have (random)
% variable delays up to a given value for that simulation.
%
% We also performed a set of experiments that are cated for larger
% deployments, similar to ones present in modern data centers. These assume
% that the network size can be thousands of nodes (we test up to 10k), have
% low delays (we test up to 5), and the graph diameter is assumed to be
% small (i.e.: less than 5). To run these experiments, please see the
% comments on how to tweak the test parameters in order to execute it.
% Please beware that it requires a considerable amount of RAM to run; 
% concretely, in my machine it required more than 63GB of available 
% memory to successfully complete.
%
% For more details, please see either the README.md or our paper which can
% be found here: https://arxiv.org/abs/2101.06139
%
% Authors:
%
%  - Andreas A. Grammenos (ag926@cl.cam.ac.uk)
%  - Themistoklis Charalambous (themistoklis.charalambous@aalto.fi)
%
% License: GPLv3
%

%% Initialisation
%
clear; clc; close all;

% for reproducibility, fix the seed
rng("default")

%% Configuration
%

% nodes for regular testing
nodes_to_test = [20, 50, 100, 150, 200, 400, 600];
% nodes for "large scale" testing
% nodes_to_test = [20, 200, 500, 1000, 5000, 10000];
node_test_len = length(nodes_to_test);
% set the node length as a convenience variable
node_len_array = 1:node_test_len;

% the number of maximum iterations
max_iter = 4000;

% the number of trials to run
% trials for regular testing
trials = 10;
% trials for large scale testing
% trials = 5;
trials_arr = 1:trials;

% maximum delay bound to test
% regular testing
delay_to_test = [1, 5, 10, 15, 20, 30];
% large scale testing
% delay_to_test = [1, 2, 3, 4, 5];
delay_len = length(delay_to_test);
delay_len_array = 1:length(delay_to_test);               

% the minimum workload value
workload_min = 1;
% the maximum workload value
workload_max = 2;

% check if we use variable capacities in the nodes
use_variable_capacities = 0;

% converge threshold - should we increase?
epsilon = 1e-5;

% timers
total_time_global = zeros(node_test_len, delay_len, trials);

% multiplier from sec to ns & μs
ns_mul = 1e+9;
micro_mul = 1e+6;

% converge statistics
cov_global = zeros(node_test_len, delay_len, trials);

cov_flips_global = zeros(node_test_len, delay_len, trials);

% converge statistics
cov_min_global = zeros(node_test_len, delay_len, trials);
cov_max_global = zeros(node_test_len, delay_len, trials);
cov_mean_global = zeros(node_test_len, delay_len, trials);
cov_win_global = zeros(node_test_len, delay_len, trials);

% setup variables
params.type = "async";          % normal async
params.type = "async-dc";       % data center async
params.pflag = 1;               % enable printing
params = setup_vars(params);    % setup environment variables

% which figures to print (can be quite spammy)
print_z_ratios = 1;             % print the z-ratios
print_all_ratios = 0;           % print all z, x, and y graphs

%% Run the trials
%

% run for the specific number of trials
for t=trials_arr
  fprintf("\n** Running trial %d\n", t);
  
  for d=delay_len_array
    
    % current maximum delay
    max_delay = delay_to_test(d);
    fprintf(" == Running for delay %d (trial %d/%d)\n", max_delay, t, trials); 
    
    % run for the nodes
    for n=node_len_array
      fprintf("  ++ Running for %d nodes (delay: %d, trial %d/%d)\n", ...
        nodes_to_test(n), max_delay, t, trials); 
      % generate the graph
      [~, diameter, nodes, AdjMatrix] = gen_graph(nodes_to_test(n));
      
      % get the node range
      graph_nodes_array = 1:nodes;
      
      outdegree_plus_one = sum(AdjMatrix);
      P = AdjMatrix * diag(1 ./ outdegree_plus_one);

      % generate load vector based on the min/max workload range.
      x_0 = gen_workload(workload_min, workload_max, nodes);
      
      % create the node capacities
      if use_variable_capacities == 0
        y_0 = ones(nodes, 1);
      else
        y1= 3 * ones(ceil(nodes / 2), 1);
        y2= 5 * ones(floor(nodes / 2), 1);
        y_0 = [y1; y2];
      end

      z_0 = x_0./y_0;
      
      xd_k = [x_0' zeros(1, max_delay * nodes)]'; xd_arxiv = x_0;
      yd_k = [y_0' zeros(1, max_delay * nodes)]'; yd_arxiv = y_0;
      zd_k = z_0; zd_arxiv = z_0;

      can_terminate = 0;
      average_X = mean(x_0);
      flag = zeros(nodes, 1);

      M = max(x_0) * ones(nodes, 1);
      m = min(x_0) * ones(nodes, 1);

      M_rxiv = [0 max(x_0)];
      m_rxiv = [0 min(x_0)];

      % the converged stats, the structure is as follows:
      %
      % row count: nodes
      %
      % description of each indice:
      %
      % 1. min iteration that satisfied the constraint
      % 2. max iteration that satisfied the contraint (and then not
      % changed)
      % 3. number of flips between min, max
      %
      %
      node_stats = zeros(nodes, 3);
      nodes_converged = 0;
      
      % reset the counter
      i = 1;
      
      % compute the mod factor for checking
      mod_factor = (1 + max_delay) * diameter;
      
      % start ticking the global one
      outer_tic = tic;
      gen_tic_cnt = 0;
      % run the consensus loop
      while i < max_iter && can_terminate == 0 

        % -- consensus inner loop -- %
        
        for j=graph_nodes_array
          if flag(j) == 0
            % make sure its the correct time to check
            if mod(i, mod_factor) == 0 && i > 1
              % check if the nodes are within the converge threshold for
              % the round; if so, raise the flag for that node
              if M(j) - m(j) < epsilon
                  prev_flag = flag(j);
                  flag(j) = 1;
                  fprintf("\t -- DEBUG: Node %d satisfied contraint at time step %d.\n", j, i)
                  
                  % -- adjust node stats -- %
                  
                  % if min is not set and the node just converged, set it.
                  if node_stats(j, 1) == 0
                    node_stats(j, 1) = i;
                  end
                  
                  % set the max only if was flipped - it is also
                  % monotonically increasing hence at every check i will be
                  % greater.
                  if prev_flag == 0
                    % change the node max timestep
                    node_stats(j, 2) = i;
                    % this counts as a flip, so advance it.
                    node_stats(j, 3) = node_stats(j, 3) + 1;
                  end
                  
                  % -- end adjust node status -- %
                  
              end
              
              M(j) = bar_zd_k(j);
              m(j) = bar_zd_k(j);
            end
          else
            % in case of M(j) - m(j) >= e, while previously satisfied: flip!
            if M(j) - m(j) >= epsilon
              fprintf("\t -- DEBUG: Node %d at iter %d does not satisfy the constraint.\n", j, i)
              % flip the node - actual flips are counted in in the top loop
              flag(j) = 0;             
            end
            % in case of M(j) - m(j) < e, do nothing
          end
          % fprintf("i: %d, j: %d, M(j): %d, m(j): %d\n", i, j, M(j), m(j));
          if ~isnan(M(j)) && ~isnan(m(j))
            assert(M(j) >= m(j) || (abs(M(j) - m(j)) < epsilon));
          end   
          
        end
        
        % -- end consensus inner loop -- %

        % -- state matrix start -- %

        % Create matrix A
        gen_tic = tic;
        [A] = gen_state_matrix(P, max_delay);
        gen_toc = toc(gen_tic);
        gen_tic_cnt = gen_tic_cnt + gen_toc;

        % Iterations with delays for X
        xd_k = A * xd_k; 
        bar_xd_k = xd_k(graph_nodes_array); 
        xd_arxiv = [xd_arxiv bar_xd_k];

        % Iterations with delays for Y
        yd_k = A * yd_k; 
        bar_yd_k = yd_k(graph_nodes_array);
        yd_arxiv = [yd_arxiv bar_yd_k];

        % Ratio of values: Z=X/Y.
        bar_zd_k = bar_xd_k ./ bar_yd_k;
        zd_arxiv = [zd_arxiv bar_zd_k];
        
        % change all nan to zero, if any - for numerical stability
        % bar_zd_k(isnan(bar_zd_k)) = 0;

        % check if a 2 x diameter of the graph plus avg delay has passed in
        % terms of iterations
        if mod(i, mod_factor) == 0 % swapped 2 * diameter * D
            M_rxiv = [M_rxiv; i max(bar_zd_k)];
            m_rxiv = [m_rxiv; i min(bar_zd_k)];
        end

        % find the minimum and maximum of B and store it.
        for j=graph_nodes_array
            B = AdjMatrix(j,:) .* bar_zd_k';
            M(j) = max(B);
            b_min = B(B ~= 0);
            % handle edge case
            if isempty(b_min)
              m(j) = 0;
            else
              m(j) = min(b_min);
            end
        end

        % -- state matrix end -- %

        % check the termination condition, which is satisfied once all 
        % nodes have raised the termination flag
        if sum(flag) == nodes
          can_terminate = 1;
          nodes_converged = i;
        else
          i = i + 1;
        end
        % - end while
      end
       
      % stop ticking the global loop
      total_time_global(n, d, t) = toc(outer_tic) - gen_tic_cnt;
      
      
      % -- compute aggregate statistics for the parameter tuple -- %
      
      % compute the stats
      cov_min = min(node_stats(:, 1));        % min converge.
      cov_max = max(node_stats(:, 2));        % max converge.
      cov_mean = mean(node_stats(:, 2));      % mean converge for max
      flip_per = sum(node_stats(:, 3))/nodes; % flip percentage.
      
      % update the variables.
      cov_min_global(n, d, t) = cov_min;
      cov_max_global(n, d, t) = cov_max;
      cov_mean_global(n, d, t) = cov_mean;
      cov_win_global(n, d, t) = cov_max - cov_min; 
      assert(cov_win_global(n, d, t) >= 0)
      cov_flips_global(n, d, t) = cov_flips_global(n, d, t) + flip_per;
      
      % -- end aggregate statistics for the parameter tuple -- %
      
      % print one sample - if enabled
      if t == 1 && print_z_ratios == 1
        
        % get the size of the zd
        [x_len, y_len] = size(zd_arxiv);
        % adjust the iteration len
        iter_len = y_len - 1;
        
        % fix a corner case where rows are equal to columns
        if x_len == x_len
          xd_arxiv = xd_arxiv';
          yd_arxiv = yd_arxiv';
          zd_arxiv = zd_arxiv';
        end
        
        % print both the x,y ratios as well.
        if print_all_ratios == 1
          fig = figure;
          hold on; box on;
            plot(0:iter_len, xd_arxiv);
            axis([0 iter_len min(xd_arxiv, [], 'all')-.2 max(xd_arxiv, [], 'all')+.2])
            set(gca, "FontSize", 20, "FontName", "Times New Roman");
            xl = xlabel("Iterations", "FontSize", 20, "FontName", "Times New Roman");
            set(xl, "Interpreter", "Latex");
            yl=ylabel('$X$ Values', "FontSize", 20, "FontName", "Times New Roman");
            set(yl, "Interpreter", "Latex");
            title(sprintf("Nodes %d", nodes))
          hold off;
          % print the figure 
          st = sprintf("nmax_%d_trials_%d_achieve_bar_x", nodes_to_test(end), trials);
          print_fig(fig, st, params);

          fig = figure;
          hold on; box on;
            plot(0:iter_len, yd_arxiv);
            axis([0 iter_len min(yd_arxiv, [], 'all')-.2 max(yd_arxiv, [], 'all')+.2])
            set(gca, "FontSize", 20, "FontName", "Times New Roman");
            xl = xlabel("Iterations", "FontSize", 20, "FontName", "Times New Roman");
            set(xl, "Interpreter", "Latex");
            yl=ylabel('$Y$ Values', "FontSize", 20, "FontName", "Times New Roman");
            set(yl, "Interpreter", "Latex");
            title(sprintf("Nodes %d", nodes))
          hold off;
          % print the figure 
          st = sprintf("nmax_%d_trials_%d_achieve_bar_y", nodes_to_test(end), trials);
          print_fig(fig, st, params);
        end

        fig = figure; 
        hold on; box on;
          plot(0:iter_len, zd_arxiv);
          axis([0 iter_len min(zd_arxiv, [], 'all')-.2 max(zd_arxiv, [], 'all')+.2])
          set(gca, "FontSize", 20, "FontName", "Times New Roman");
          xl = xlabel("Iterations", "FontSize", 20, "FontName", "Times New Roman");
          set(xl,"Interpreter","Latex");
          yl = ylabel("Ratios", "FontSize", 20, "FontName", "Times New Roman");
          set(yl, "Interpreter", "Latex");
          plot([0, iter_len], [average_X, average_X], 'k--')
          stairs(M_rxiv(:, 1), M_rxiv(:, 2), 'r--')
          stairs(m_rxiv(:, 1), m_rxiv(:, 2), 'b--')
          title(sprintf("Nodes %d, Delay %d, Trial: %d, Diameter: %d", ...
            nodes, max_delay, t, diameter))
        hold off;
        % print the figure 
        st = sprintf("nmax_%d_delay_%d_trials_%d_achieve_bar_z", ...
          nodes, max_delay, t);
        print_fig(fig, st, params);
      end
      
      % check if we terminated or exhausted the max iterations
      if can_terminate == 0
        fprintf("   ## WARN: Termination condition was NOT met using max iter of %d\n", i);
      else
        fprintf("   ^^ INFO: All nodes converged at iteration %d.\n", i)
        cov_global(n, d, t) = i;
      end
      fprintf("  ++ Finished running for %d nodes (delay: %d, trial %d/%d)\n", ...
        nodes_to_test(n), max_delay, t, trials); 
    end
    
    fprintf(" == Finished running for delay %d (trial %d/%d)\n", max_delay, t, trials);
  end
  fprintf("** Finished trial %d\n", t);
end

%% Plot aggregate figures
%

% -- plot delay -- %

fig = figure;
hold on; box on;
wlegs = cell(1, node_test_len);
for i=node_len_array
    cur = squeeze(cov_global(i, :, :));
    err_std = std(cur, 0, 2)' / sqrt(trials);
    errorbar(delay_len_array, mean(cur, 2)', err_std, "LineWidth", 2);
    wlegs{i} = sprintf("N=%d", nodes_to_test(i));
end
hold off;

ylabel("Iterations", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xlabel("Delay ($\tau$)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
title("Converge iterations per delay used", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(delay_len_array)
xticklabels(num2cell(delay_to_test))
legend(wlegs);

% print the figure 
st = sprintf("nmax_%d_trials_%d_delay_scaling", nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- plot flip percentage -- %
fig = figure;
hold on; box on;
wlegs = cell(1, node_test_len);
for i=node_len_array
    cur = squeeze(cov_flips_global(i, :, :));
    err_std = std(cur, 0, 2)' / sqrt(trials);
    errorbar(delay_len_array, mean(cur, 2)', err_std, "LineWidth", 2);
    wlegs{i} = sprintf("N=%d", nodes_to_test(i));
end
hold off;

ylabel("Flips \% over 1", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xlabel("Delay ($\tau$)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
title("Flip percentage during iterations convergence", ... 
  "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(delay_len_array)
xticklabels(num2cell(delay_to_test))
legend(wlegs);

% print the figure 
st = sprintf("nmax_%d_trials_%d_flip_percentage", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- plot flip statistics -- %

wlegs = ["min", "max", "mean", "window"];
msg = "Converge statistics for Delay ($\\tau$): %d";
for j=delay_len_array
  % create the figure
  fig = figure;
  
  % enable hold and boxing
  hold on; box on;
  
  % min cov
  cur = squeeze(cov_min_global(:, j, :));
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % max cov
  cur = squeeze(cov_max_global(:, j, :));
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % mean
  cur = squeeze(cov_mean_global(:, j, :));
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % window cov
  cur = squeeze(cov_win_global(:, j, :));
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  hold off;
  
  % handle figure presentation
  ylabel("Iterations", "FontName", "Times New Roman");
  xlabel("Nodes", "Interpreter", "Latex", "FontName", "Times New Roman");
  st = sprintf(msg, delay_to_test(j));
  title(st, "Interpreter", "Latex", "FontName", "Times New Roman");
  xticks(node_len_array)
  xticklabels(num2cell(nodes_to_test))
  legend(wlegs);
  
  % print the figure 
  st = sprintf("nmax_%d_trials_%d_delay_%d_converge_statistics", ...
    nodes_to_test(end), trials, delay_to_test(j));
  print_fig(fig, st, params);
end

% -- plot execution time -- %

fig = figure;
hold on; box on;
wlegs = cell(1, node_test_len);
for i=node_len_array
    cur = squeeze(total_time_global(i, :, :));
    err_std = std(cur, 0, 2)' / sqrt(trials);
    errorbar(delay_len_array, mean(cur, 2)', err_std, "LineWidth", 2);
    wlegs{i} = sprintf("N=%d", nodes_to_test(i));
end
hold off;

ylabel("Time (s)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xlabel("Delay ($\tau$)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
title("Total time to converge", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(delay_len_array)
xticklabels(num2cell(delay_to_test))
legend(wlegs);

% print the figure 
st = sprintf("nmax_%d_trials_%d_total_exec_time", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);