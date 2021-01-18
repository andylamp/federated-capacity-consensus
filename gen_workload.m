function [workload] = gen_workload(min_w, max_w, len)
%GEN_WORKLOAD Generates the workload vector based on min/max values
%
  % check if the minimum workload is greater than the maximum
  if min_w > max_w
    error("Minimum workload cannot be greater than maximum")
  end
  
  % check if the minimum workload is less than zero
  if min_w < 0
    error("Minimum workload cannot be zero")
  end
  
  % return the workload vector in row-major format
  workload = randi([min_w, max_w], len, 1);
end

