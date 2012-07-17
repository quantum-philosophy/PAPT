function [ nodes, weights, count ] = doConstructGrid(sdim, order)
  %
  % Description:
  %
  %   Construct a sparse grid based on the Smolyak algorithm with
  %   the 1D Gauss-Hermite rule inside.
  %
  %   NOTE: The physicists' type of the Hermite polynomials is assumed,
  %   hence, the weight function is exp(-x^2).
  %

  %
  % We have integrals in the form <f(Z), Psi_i(Z)> where the maximal order
  % of Psi_i(Z) is `order'. At the same time, we are trying to approximate
  % f(Z) with a polynomial of the same order `order'. Therefore, the inner
  % product can be a integral of a (2 * `order')-order polynomial.
  %
  % A Gaussian quadrature rule with n points is exact of (2 * n - 1)-order
  % polynomials, hence, we need to take (n + 1) points to cover
  % (2 * `order')-order polynomials as well.
  %
  % So, with the tensor product everything is clear.
  %
  TPcount1D = order + 1;
  TPcount = TPcount1D^sdim;

  %
  % For sparse grids, however, there are no such guarantees^. Instead,
  % people are talking in terms of `levels' (of accuracy), which are
  % simply indexes of families of sparse grids [1].
  %
  % [1] http://people.sc.fsu.edu/~jburkardt/m_src/sparse_grid_hermite/sparse_grid_hermite.html
  %
  % Let this level be the order of the PC expansion.
  %
  SGlevel = order;
  SGcount = sparse_grid_herm_size(sdim, SGlevel);

  count = min(TPcount, SGcount);

  warn(sdim, order, count, TPcount < SGcount);

  if TPcount < SGcount
    [ nodes1d, weights1d ] = hermite_compute(TPcount1D);

    nodes = {};
    weights = {};

    for i = 1:sdim
      nodes{end + 1} = nodes1d;
      weights{end + 1} = weights1d;
    end

    [ nodes, weights ] = tensor_product(nodes, weights);

    nodes = transpose(nodes);
    weights = transpose(weights);
  else
    [ weights, nodes ] = sparse_grid_herm(sdim, SGlevel, SGcount);
  end
end

function warn(sdim, order, count, TP)
  if TP
    method = 'Tensor product';
  else
    method = 'Smolyak';
  end

  debug('------------------------------\n');
  debug('Constructing a new grid:\n');
  debug('  Stochastic dimension: %d\n', sdim);
  debug('  Polynomial order: %d\n', order);
  debug('  Number of points: %d\n', count);
  debug('  Method: %s\n', method);
  debug('------------------------------\n');
end
