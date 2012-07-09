function [ x, psi, norm, count ] = prepareExpansion(dimension, order, gq)
  %
  % Description:
  %
  %   Calculate multivariate (Hermite) polynomials of `dimension' r.v.'s
  %   of the maximal order `order'. In addition, compute the variances
  %   of each of the polynomials.
  %
  % Input:
  %
  %   * dimension - the stochastic dimension of the expansion.
  %   * order     - the maximal order of the polynomials.
  %   * gq        - the Gauss quadrature to integrate with.
  %
  % Output:
  %
  %   * x     - a vector of `dimension' symbolic variables.
  %   * psi   - a vector of the Hermite polynomials.
  %   * norm  - a vector of the normalization coefficients of the expansion.
  %   * count - the total number of the polynomials in `psi'.
  %

  filename = [ 'PC_d', num2str(dimension), '_o', num2str(order), '.mat' ];
  filename = Utils.resolvePath(filename);

  if exist(filename, 'file')
    load(filename);
  else
    count = PC.calculateCount(dimension, order);

    for i = 1:dimension
      x(i) = sym([ 'x', num2str(i) ], 'real');
    end

    psi = PC.constructXD(x, count);

    norm = zeros(1, count);
    norm(1) = 1;
    for i = 2:count
      norm(i) = gq.integrate(@(s) subs(psi(i) * psi(i), x, s));
    end

    save(filename, 'x', 'psi', 'norm', 'count');
  end
end
