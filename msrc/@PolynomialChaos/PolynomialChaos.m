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
    % The Gaussian quadrature (Gauss-Hermite) for the purpose of
    % numerical integration. It is an optimized version, an instance
    % of GaussianQuadrature.Chaos.
    %
    qd

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
    function pc = PolynomialChaos(sdim, order)
      if nargin < 2, order = 2; end

      pc.sdim = sdim;
      pc.order = order;

      [ pc.x, pc.psi, index ] = pc.prepareExpansion(sdim, order);

      pc.qd = GaussianQuadrature.MultiProbabilists(pc.x, pc.psi, index);

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

  methods (Static)
    function count = countTerms(sdim, order, rule)
      if nargin < 3, rule = 'TD'; end

      switch (upper(rule))
      case 'TD'
        count = factorial(sdim + order) / factorial(sdim) / factorial(order);
      case 'TP'
        count = (1 + order)^sdim;
      otherwise
        error('Unknown rule.');
      end
    end
  end

  methods (Static, Access = 'private')
    psi = construct1D(x, order);
    [ psi, index ] = constructMD(x, order);
    [ x, psi, index ] = doPrepareExpansion(sdim, order);

    function [ x, psi, index ] = prepareExpansion(sdim, order)
      %
      % A wrapper to cache the result of `doPrepareExpansion'.
      %

      filename = [ 'PolynomialChaos_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ x, psi, index ] = PolynomialChaos.doPrepareExpansion(sdim, order);
        save(filename, 'x', 'psi', 'index');
      end
    end
  end
end
