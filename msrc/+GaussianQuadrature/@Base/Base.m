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
    function gq = Base(x, psi, order);
      [ gq.nodes, gq.plainGrid, gq.niceGrid, gq.norm ] = gq.precomputeGrid(x, psi, order);
      gq.points = size(gq.nodes, 2);
    end
  end

  methods (Static, Access = 'protected')
    function [ nodes, plainGrid, niceGrid, norm ] = precomputeGrid(x, psi, order)
      error('To be implemented.');
    end
  end
end
