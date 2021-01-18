function [G, distance, nodes, adj] = gen_graph(nodes, opt)
%GEN_GRAPH helper function that generates our graph
%

  % let the user know of the mishap.
  if nodes < 1
    error("we requite at least 1 node")
  end
  
  % if opt was not supplied initialise an empty struct
  if nargin < 2
    opt = struct;
  end
  
  % check if we want to use a digraph or not
  if ~isfield(opt, "is_digraph")
    opt.is_digraph = 1;
  end
  
  % do not show the graphs by default
  if ~isfield(opt, "show_graph")
    opt.show_graph = 0;
  end
  
  % check if we are provided with an adjacency matrix
  if isfield(opt, "adj")
    opt.adj_matrix_method = "provided";
  elseif ~isfield(opt, "adj_matrix_method")
    % if not provided, use the "fancy method"
    opt.adj_matrix_method = "fancy";
    % opt.adj_matrix_method = "regular";
  end
  
  % check if we have a connectivity constraint for the generated graph
  if ~isfield(opt, "connectivity")
    % enforce strong connectivity requirement, if we have a digraph and 
    % was not provided
    if opt.is_digraph == 1
      opt.connectivity = "strong";
    else
      opt.connectivity = "dontcare";
    end
  end
  
  % check if we have a max graph generation iteration constraint
  if ~isfield(opt, "max_graph_iter")
    opt.max_graph_iter = 100000;
  end
  
  % check if we require a min diameter for our graph
  if ~isfield(opt, "min_graph_diameter")
    opt.min_graph_diameter = 1;
  end
 
  has_min_diameter = 0;
  cnt = 0;
  while cnt < opt.max_graph_iter
    % generate a random binary graph adjacency matrix, 
    % which is also returned
    if opt.adj_matrix_method == "fancy"
      % specify some exact number of zeros
      zero_number = ceil(0.85 * (nodes^2));
      % create the matrix using full zeros
      adj = ones(nodes);
      % now set the number of zeros
      adj(1:zero_number) = 0; 
      adj(randperm(numel(adj))) = adj;
    elseif opt.adj_matrix_method == "regular"
      % a more basic method
      adj = round(rand(nodes));
      adj = triu(adj) + triu(adj, 1)';
      adj = adj - diag(diag(adj));
    end

    % now create a graph g, based on either the provided or generated 
    % adjacency matrix
    if opt.is_digraph == 1
      G = digraph(adj);
    else
      G = graph(adj, 'upper');
    end
    
    % ensure we have the proper diameter, first compute its distance
    distance = max(max(distances(G), [], 'omitnan'));
    
    % now raise the flag that we have the minimum required distance
    if ~isinf(distance) && distance >= opt.min_graph_diameter
      has_min_diameter = 1;
    end

    % ensure we have the proper connectivity
    if opt.connectivity == "strong"
      bins = conncomp(G, 'Type', 'strong');
      connecticity_satisfied = all(bins == 1);
    elseif opt.connectivity == "weak"
      bins = conncomp(G, 'Type', 'weak');
      connecticity_satisfied = all(bins == 1);
    elseif opt.connectivity == "dontcare"
      % do nothing
      connecticity_satisfied = 1;
    else
      error("Invalid connectivity type provided.")
    end
    
    % check if both connectivity and min distance is satisfied
    if connecticity_satisfied == 1 && has_min_diameter == 1
      break
    end
    
    % increase the counter
    cnt = cnt + 1;
  end
  
  if cnt == opt.max_graph_iter
    if has_min_diameter == 0
      fprintf("min graph distance was not satisfied\n");
    end
    error("could not converge for %d iterations", opt.max_graph_iter);
  end
  
  fprintf("\t$$ Graph gen took %d iterations, diameter %d - connection type: %s\n", ...
    cnt, distance, opt.connectivity);
  
  % check if we show the graph
  if opt.show_graph == 1
    figure;
    plot(G);
    title(sprintf("Graph with %d nodes", nodes));
  end

end

