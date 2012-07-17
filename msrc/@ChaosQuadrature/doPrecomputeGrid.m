function [ nodes, plainGrid, niceGrid ] = doPrecomputeGrid(x, psi, order)
  sdim = length(x);
  terms = length(psi);

  [ nodes, weights, count ] = GaussianQuadrature.constructGrid(sdim, order);

  warn(sdim, order, terms, count);

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

function warn(sdim, order, terms, count)
  debug('------------------------------\n');
  debug('Precomputing a new PC grid:\n');
  debug('  Stochastic dimension: %d\n', sdim);
  debug('  Polynomial order: %d\n', order);
  debug('  Number of terms: %d\n', terms);
  debug('  Number of points: %d\n', count);
  debug('------------------------------\n');
end
