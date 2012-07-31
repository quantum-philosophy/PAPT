classdef Base < handle
  properties (SetAccess = 'private')
    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid.
    %
    plainGrid

    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid multiplied by the
    % corresponding weight.
    %
    niceGrid

    %
    % The normalization coefficients of the expansion, i.e., <psi_i^2>.
    %
    norm
  end

  properties (SetAccess = 'private')
    %
    % The evaluation points for the integration.
    %
    nodes

    %
    % The number of points.
    %
    points
  end

  methods
    function gq = Base(x, psi, order, varargin)
      [ gq.nodes, gq.plainGrid, gq.niceGrid, gq.norm ] = gq.precomputeGrid(x, psi, order, varargin{:});
      gq.points = size(gq.nodes, 2);
    end
  end

  methods (Access = 'protected')
    [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(gq, x, psi, order, varargin)
  end

  methods (Access = 'protected')
    function quadratureOrder = polynomialOrderToQuadratureOrder(gq, polynomialOrder)
      %
      % The order of a Gaussian quadrature rule, denoted by `order', means that
      % the rule is exact for (2 * `order' - 1)-order polynomials. We want to have
      % exactness for polynomials of order (2 * `order'), therefore, the rule order
      % should be increased by one.
      %
      quadratureOrder = polynomialOrder + 1;
    end

    function points = countTensorProductPoints(gq, sdim, order)
      points = order^sdim;
    end

    function [ nodes, plainGrid, niceGrid, norm ] = precomputeGrid(gq, x, psi, order, varargin)
      sdim = length(x);

      filename = [ class(gq), '_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, plainGrid, niceGrid, norm ] = gq.doPrecomputeGrid(x, psi, order, varargin{:});
        save(filename, 'nodes', 'plainGrid', 'niceGrid', 'norm');
      end
    end
  end
end
