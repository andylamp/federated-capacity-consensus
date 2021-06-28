function [out_mat] = interp_delay(in_mat)
%INTERP_DELAY Interpolate values in case of nan's due to delay.
% This is performed primarity to have a nice printed output.
  % check if the input vector has nan
  nans = isnan(in_mat);
  
  % check if we have all values, if we do then just return.
  if ~any(nans)
    out_mat = in_mat;
    return
  end
  
  % use the average of the previous and next values to interpolate
  out_mat = .5 * (fillmissing(in_mat, 'previous') + ...
    fillmissing(in_mat, 'next'));
end

