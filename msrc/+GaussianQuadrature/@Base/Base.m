classdef Base < handle
  %
  % (*) The quadrature rule that is used here (Gauss-Hermite) assumes
  % the physicists' weight function. We need the probabilists' one,
  % therefore, a proper conversion is to be performed.
  %
  properties (SetAccess = 'private')
    %
    % The dimension of the grid, i.e., the number of variables.
    %
    sdim

    %
    % The order of the polynomial expansion.
    %
    order

    %
    % The total number of nodes.
    %
    count

    %
    % The nodes where integrands will be evaluated.
    %
    nodes

    %
    % The corresponding weights.
    %
    weights
  end

  methods
    function gq = Base(sdim, order)
      gq.sdim = sdim;
      gq.order = order;

      [ nodes, gq.weights, gq.count ] = ...
        GaussianQuadrature.constructGrid(sdim, order);

      %
      % See (*) to justify the need of sqrt(2).
      %
      gq.nodes = nodes * sqrt(2);
    end

    function result = integrate(gq, f, ddim)
      count = gq.count;
      nodes = gq.nodes;

      samples = zeros(ddim, count);

      for i = 1:count
        samples(:, i) = f(nodes(:, i));
      end

      %
      % See (*) to understand why we do not need 2 next to pi here.
      %
      result = sum(samples .* irep(gq.weights, ddim, 1), 2) ./ pi^(gq.sdim / 2);
    end
  end
end
