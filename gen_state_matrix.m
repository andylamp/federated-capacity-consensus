function [P_aug] = gen_state_matrix(P, max_delay)
%% Creates the state matrix A
%

  % get the size
  [dim, N] = size(P);
 
  % check for the validity
  if dim ~= N
    error("This function is designed to work with square matrices")
  end
  
  % construct the padding matrices
  I = eye(dim);   % the identity pad
  O = zeros(dim); % the zero pad
  
  % generate the random states for all the different delays
  R = max_delay * rand(dim); 
  % copy matrix
  P0 = P; 
  % now filter all non-zero values and set them to 1
  P0(P ~= 0) = 1;
  % 
  P0= P0 - diag(diag(P0));
  
  R = R .* P0;
  % round everything up to an integer
  R = round(R);
  
  % generate the tables
  P_i = zeros(max_delay + 1, dim, dim);
  
  for i=1:dim
      for j=1:dim
        % loop for each possible delay
        for d=1:max_delay+1
          % check if we matched our delay
          if R(i, j) == d-1
            % assign the value to to the state based on P
            P_i(d, i, j) = P(i, j);
            continue
          end
        end
      end
  end

  % now construct the state matrix based on the above - it is much faster
  % as it uses preallocation and direct indices to perform the mapping.
  P_aug = zeros(dim * (max_delay + 1));
  % construct rows
  for i=1:max_delay+1
    % assign columns
    i_start = (i-1)*dim + 1;
    i_end = i*dim;
    for j=1:max_delay+1
      j_start = (j-1)*dim + 1;
      j_end = j*dim;
      if j == 1
        P_aug(i_start:i_end, j_start:j_end) = squeeze(P_i(i, :, :)); 
      elseif j == i+1
        P_aug(i_start:i_end, j_start:j_end) = I;
      else
        P_aug(i_start:i_end, j_start:j_end) = O;
      end
    end
  end
end

