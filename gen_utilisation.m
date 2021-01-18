function [utilisation] = gen_utilisation(nodes, use_round)
%GEN_UTILISATION Generates the booking utilisation vector
%
  % check if round was supplied
  if nargin < 2
    use_round = 1;
  end
  
  % generate the utilisation vector
  utilisation = rand(nodes, 1);
  
  % if we use round, then we round to first decimal - i.e. 0.12 -> 0.1
  if use_round == 1
    utilisation = round(utilisation, 1);
  end
end

