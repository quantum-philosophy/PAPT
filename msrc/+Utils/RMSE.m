function rmse = RMSE(observed, predicted)
  o = observed(:);
  p = predicted(:);
  nrmse = sqrt(sum((o - p) .^ 2) / numel(o));
end
