function [ nodes, plainGrid, niceGrid ] = doPrecomputeGrid(x, psi, level)

  sdim = length(x);
  terms = length(psi);

  [ nodes, weights, count ] = GaussianQuadrature.constructGrid(sdim, level);

  %
  % See (*) in GaussianQuadrature.
  %
  nodes = nodes * sqrt(2);

  plainGrid = zeros(terms, count);
  niceGrid = zeros(terms, count);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');

    plainGrid(i, :) = f(nodes);

    %
    % See (*) in GaussianQuadrature.
    %
    niceGrid(i, :) = plainGrid(i, :) .* weights ./ pi^(sdim / 2);
  end
end
