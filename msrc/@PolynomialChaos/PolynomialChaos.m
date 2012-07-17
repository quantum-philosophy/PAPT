classdef PolynomialChaos < handle
  properties (SetAccess = 'protected')
    %
    % The stochastic dimension of the polynomial chaos (PC) expansion,
    % i.e., the number of r.v.'s involved.
    %
    sdim

    %
    % The maximal order of the multivariate (Hermite) polynomials.
    %
    order

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
    % The Gaussian quadrature (Gauss-Hermite) for the purpose of
    % numerical integration. It is an optimized version, an instance
    % of ChaosGuadrature.
    %
    cq

    %
    % The total number of polynomials in the expansion.
    %
    terms

    %
    % The function that takes two vector-valued arguments and
    % evaluates the whole PC expansion.
    %
    evaluator
  end

  methods
    function pc = PolynomialChaos(sdim, order, level)
      if nargin < 2, order = 2; end
      if nargin < 3, level = ceil((order + 1)/2); end

      pc.sdim = sdim;
      pc.order = order;

      [ pc.x, pc.psi, pc.norm, pc.cq ] = ...
        pc.prepareExpansion(sdim, order, level);

      pc.terms = length(pc.psi);

      %
      % Construct the final polynomial with abstract coefficients
      % and produce a function to evaluate it.
      %
      for i = 1:pc.terms
        a(i) = ipoly([ 'a', num2str(i) ]);
      end
      pc.evaluator = Utils.toFunction(sum(a .* pc.psi), pc.x, 'rows', a);
    end
  end

  methods (Static, Access = 'private')
    psi = construct1D(x, terms);
    psi = constructXD(x, terms);
    [ x, psi, norm, cq ] = doPrepareExpansion(sdim, order, level);

    function terms = countTerms(sdim, order)
      terms = factorial(sdim + order) / (factorial(sdim) * factorial(order));
    end

    function [ x, psi, norm, cq ] = prepareExpansion(sdim, order, level)
      %
      % A wrapper to cache the result of `doPrepareExpansion'.
      %
      % NOTE: We do not cache `cq' explicitly.
      %

      filename = [ 'PC_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
        cq = ChaosQuadrature(x, psi, order, level);
      else
        [ x, psi, norm, cq ] = ...
          PolynomialChaos.doPrepareExpansion(sdim, order, level);
        save(filename, 'x', 'psi', 'norm');
      end
    end
  end
end
