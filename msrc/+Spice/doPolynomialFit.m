function [ fitresult, gof, LTmean, LTstd ] = doFit(L, T, I, order)
  if nargin < 4, order = [ 2 2 ]; end

  [ X, Y, Z ] = prepareSurfaceData(L, T, I);
  X = [ X, Y ];
  Y = Z;

  ft = fittype(sprintf('poly%d%d', order(1), order(2)));
  opts = fitoptions(ft);

  opts.Normalize = 'off';

  [ X, LTmean, LTstd ] = curvefit.normalize(X);

  vars = order(1) * order(2) + 1;
  opts.Lower = ones(1, vars) * -Inf;
  opts.Upper = ones(1, vars) * Inf;

  [ fitresult, gof ] = fit(X, Y, ft, opts);
end
