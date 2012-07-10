classdef PolynomialChaos < handle
  properties (SetAccess = 'protected')
    %
    % The stochastic dimension of the polynomial chaos (PC) expansion,
    % i.e., the number of r.v.'s involved.
    %
    dimension

    %
    % The maximal order of the multivariate (Hermite) polynomials.
    %
    order

    %
    % The Gaussian quadrature (Gauss-Hermite) for the purpose of
    % numerical integration.
    %
    gq

    %
    % The symbolic variables used in the polynomials.
    %
    x

    %
    % The polynomials of the expansion.
    %
    psi

    %
    % The normalization coefficients of the expansion, i.e.,
    % the variances of the polynomials E<psi_i^2>.
    %
    norm

    %
    % The total number of polynomials in the expansion.
    %
    terms
  end

  methods
    function pc = PolynomialChaos(dimension, order, level)
      if nargin < 2, order = 2; end
      if nargin < 3, level = ceil((order + 1)/2); end

      pc.dimension = dimension;
      pc.order = order;

      pc.gq = GaussianQuadrature(dimension, level);

      [ pc.x, pc.psi, pc.norm, pc.terms ] = ...
        pc.prepareExpansion(dimension, order, pc.gq);
    end
  end

  methods (Static, Access = 'private')
    function terms = countTerms(dimension, order)
      terms = factorial(dimension + order) / ...
        (factorial(dimension) * factorial(order));
    end

    psi = construct1D(x, terms);
    psi = constructXD(x, terms);

    [ x, psi, norm, terms ] = prepareExpansion(dimension, order, gq);
  end
end
