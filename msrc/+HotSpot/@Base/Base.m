classdef Base < handle
  properties (SetAccess = 'private')
    %
    % The number of thermal nodes in the thermal RC circuit.
    %
    nodes

    %
    % The number of active node, i.e., those that correspond to
    % the processing elements (cores).
    %
    cores

    %
    % The sampling interval of the simulation.
    %
    dt

    %
    % The ambient temperature.
    %
    Tamb

    %
    % The stochastic dimension of the analysis, i.e., the number
    % of r.v.'s involved (obtained via the PCA).
    %
    sdim
  end

  properties (GetAccess = 'protected', SetAccess = 'private')
    %
    % The capacitance vector.
    %
    C

    %
    % The conductance matrix.
    %
    G

    %
    % The mapping matrix from the r.v.'s to the cores.
    %
    map
  end

  methods
    function hs = Base(floorplan, config, line)
      [ hs.C, hs.G, hs.nodes, hs.cores, hs.dt, hs.Tamb ] = ...
        HotSpot.Base.getCoefficients(floorplan, config, line);

      hs.nodes = size(hs.G, 1);

      %
      % Perform the PCA to extract independent r.v.'s.
      %
      [ hs.sdim, hs.map ] = PrincipalComponent.perform(hs.cores);
      assert(size(hs.map, 1) == hs.cores, 'The dimensions do not match.');
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, config, line);
  end
end