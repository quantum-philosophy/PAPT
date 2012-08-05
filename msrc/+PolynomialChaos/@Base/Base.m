classdef Base < handle
  properties (SetAccess = 'protected')
    %
    % The stochastic dimension of the polynomial chaos (PC) expansion,
    % i.e., the number of r.v.'s involved.
    %
    sdim

    %
    % The deterministic dimension of the problem.
    %
    ddim

    %
    % The maximal order of the multivariate (Hermite) polynomials.
    %
    order

    %
    % The followings comes from the Gaussian quadrature rule.
    %

    %
    % The evaluation points of the rule.
    %
    nodes

    %
    % The normalization constants of each of the polynomials, i.e., <psi_i^2>.
    %
    norm

    %
    % The grid for integration in multiple stochastic and deterministic dimensions.
    %
    grid

    %
    % The following matrices represent the expanded version of the PC
    % expansion in a sum of `mterms' monomials.
    %

    %
    % A (sdim x mterms) matrix of the exponents of the monomials.
    %
    rvPower

    %
    % A (mterms x points) matrix of precomputed values of products of r.v.'s
    % for all monomials and all points.
    %
    rvProd

    %
    % A (terms x mterms) matrix that maps the coefficients of the PC
    % expansion to the coefficients of the monomials.
    %
    coeffMap

    %
    % A shortcut for coeffMap * rvProd.
    %
    mappedRvProd

    %
    % The number of points in the rule.
    %
    points

    %
    % The total number of polynomials in the expansion.
    %
    terms
  end

  methods
    function pc = Base(dims, order, varargin)
      method = pc.prepare(varargin{:});

      pc.sdim = dims(1);
      pc.ddim = dims(2);
      pc.order = order;

      [ pc.nodes, pc.norm, pc.grid, pc.rvPower, pc.rvProd, pc.coeffMap ] = ...
        pc.prepareExpansion(pc.sdim, pc.ddim, order, method);

      pc.mappedRvProd = pc.coeffMap * pc.rvProd;

      pc.points = size(pc.nodes, 2);
      pc.terms = length(pc.norm);
    end

    function newCoeff = computeExpansion(pc, f, prevCoeff)
      grid = pc.grid;
      terms = pc.terms;

      newCoeff = zeros(pc.ddim, terms);

      if nargin > 2
        %
        % The first option is to perform a real evaluation of
        % the current PC expansion of temperature like this...
        %
        currentValue = prevCoeff * pc.mappedRvProd;
        %
        % However, it leads to an unstable behaviour of PDFs since
        % the quadratures cannot properly integrate when the leakage
        % current changes too rapidly (not smooth).
        %
        % Instead, one can try sampling the leakage current at
        % the expected temperature like this...
        %
        % currentValue = irep(prevCoeff(:, 1), 1, pc.points);
        %
        % However, in this case, expectation and variance start
        % deviating from the MCS too much. On the other hand, PDFs
        % around mean values -/+ standard deviations are pretty good.
        %
        samples = f(pc.nodes, currentValue);
      else
        samples = f(pc.nodes);
      end

      for i = 1:terms
        newCoeff(:, i) = sum(samples .* grid{i}, 2);
      end
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

    function order = indexToOrder(index)
      %
      % (-1) because in MATLAB indexes start from 1.
      %
      order = max(max(index)) - 1;
    end
  end

  methods (Abstract, Access = 'protected')
    psi = construct1D(pc, x, order);
  end

  methods (Access = 'protected')
    function nodes = generateSampleNodes(pc, points)
      nodes = normrnd(0, 1, pc.sdim, points);
    end

    function method = prepare(pc, method)
      if nargin < 2, method = struct(); end
    end
  end

  methods (Access = 'private')
    [ psi, index ] = constructMD(pc, x, order);

    function [ nodes, norm, grid, rvPower, rvProd, coeffMap ] = ...
      prepareExpansion(pc, sdim, ddim, order, method)
      %
      % A wrapper to cache the result of `doPrepareExpansion'.
      %

      filename = [ PolynomialChaos.methodStamp(method), ...
        '_', Quadrature.methodStamp(method), ...
        '_sd', num2str(sdim), ...
        '_dd', num2str(ddim), ...
        '_o', num2str(order), '.mat' ];

      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, norm, grid, rvPower, rvProd, coeffMap ] = ...
          pc.doPrepareExpansion(sdim, ddim, order, method);
        save(filename, 'nodes', 'norm', 'grid', 'rvPower', 'rvProd', 'coeffMap');
      end
    end
  end
end
