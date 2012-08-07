classdef Base < handle
  properties (Constant)
    %
    % The nominal value of the channel length and its deviation.
    %
    Lnom = 45e-9;

    %
    % The following constants are used to construct an instance
    % of the leakage model that produces `PleakPdyn' portion
    % of the given dynamic power at temperature `Tref'.
    %
    PleakPdyn = 2/3;
    Tref = Utils.toKelvin(120);
  end

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
    % The result of the PCA of the floorplan, which takes into consideration
    % the global variations as well.
    %
    pca

    %
    % The leakage model.
    %
    computeLeakage
    leakageAlpha
  end

  methods
    function hs = Base(floorplan, config, line, varargin)
      [ hs.C, hs.G, hs.nodes, hs.cores, hs.dt, hs.Tamb ] = ...
        HotSpot.Base.getCoefficients(floorplan, config, line);

      hs.nodes = size(hs.G, 1);

      %
      % Perform the PCA to extract independent r.v.'s.
      %
      [ hs.pca, hs.sdim ] = PrincipalComponent.perform(floorplan, hs.Lnom, varargin{:});
      assert(size(hs.pca, 1) == hs.cores, 'The dimensions do not match.');
    end

    function configureLeakage(hs, Pdyn)
      hs.computeLeakage = Spice.fitExponentPolynomial( ...
        'inverter_45nm', [ 1, 2 ], [ 1, 0.7, 0; 1, 1, 1 ]);

      P0 = hs.computeLeakage(hs.Lnom, hs.Tref);
      hs.leakageAlpha = hs.PleakPdyn * mean(mean(Pdyn)) / P0;
    end

    function display(hs)
      fprintf('%s:\n', class(hs));
      fprintf('  Processing elements: %d\n', hs.cores);
      fprintf('  Thermal nodes: %d\n', hs.nodes);
      fprintf('  Time step: %.2e s\n', hs.dt);
      fprintf('  Ambient temperature: %.2f C\n', Utils.toCelsius(hs.Tamb));
      fprintf('  Stochastic dimensions: %d\n', hs.sdim);
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, config, line);
  end
end
