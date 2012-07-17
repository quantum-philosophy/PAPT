classdef HotSpot < handle
  properties (SetAccess = 'protected')
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

  properties (Access = 'protected')
    %
    % Original system:
    %
    %   C * dZ/dt + G * Z = M * P
    %   T = M^T * Z + T_amb
    %
    % Transformed system:
    %
    %   dX/dt = A * X + B * P
    %   T = B^T * X + T_amb
    %

    %
    % C^(-1/2)
    %
    Cm12

    %
    % A = C^(-1/2) * (-G) * C^(-1/2)
    %
    A

    %
    % B^T = (C^(-1/2) * M)^T
    %
    BT

    %
    % Eigenvalue decomposition of Gt
    % A = Q * L * Q^T
    %
    L
    Q
    QT

    %
    % E = exp(A * t) = Q * diag(exp(li * dt)) * Q^T
    %
    E

    %
    % D = A^(-1) * (exp(A * t) - I) * B
    %   = Q * diag((exp(li * t) - 1) / li) * Q^T * B
    %
    D

    %
    % The mapping matrix from the r.v.'s to the cores.
    %
    map
  end

  methods
    function hs = HotSpot(floorplan, hsConfig, hsLine)
      [ C, G, hs.nodes, hs.cores, hs.dt, hs.Tamb ] = ...
        HotSpot.getCoefficients(floorplan, hsConfig, hsLine);

      hs.nodes = size(G, 1);

      hs.Cm12 = diag(sqrt(1 ./ C));

      M = [ diag(ones(1, hs.cores)); zeros(hs.nodes - hs.cores, hs.cores) ];

      hs.A = Utils.symmetrize(hs.Cm12 * (-G) * hs.Cm12);

      B = hs.Cm12 * M;
      hs.BT = B';

      [ Q, L ] = eig(hs.A);

      hs.L = diag(L);
      hs.Q = Q;
      hs.QT = Q';

      hs.E = hs.Q * diag(exp(hs.dt * hs.L)) * hs.QT;
      hs.D = hs.Q * diag((exp(hs.dt * hs.L) - 1) ./ hs.L) * hs.QT * B;

      %
      % Perform the PCA to extract independent r.v.'s.
      %
      [ hs.sdim, hs.map ] = PrincipalComponent.perform(hs.cores);
      assert(size(hs.map, 1) == hs.cores, 'The dimensions do not match.');
    end

    function T = solve(hs, Pdyn, rvs)
      [ cores, steps ] = size(Pdyn);
      assert(cores == hs.cores, 'The power profile is invalid.')

      %
      % General shortcuts.
      %
      nodes = hs.nodes;
      E = hs.E;
      D = hs.D;
      BT = hs.BT;
      Tamb = hs.Tamb;

      %
      % Initialize the leakage model.
      %
      sampler = Leakage(Tamb, Pdyn, hs.map);

      T = zeros(cores, steps);

      Pleak = sampler.performAtAmbient(rvs);
      X = D * (Pdyn(:, 1) + Pleak);

      for i = 2:steps
        T(:, i - 1) = BT * X + Tamb;
        Pleak = sampler.performAtGiven(T(:, i - 1), rvs);
        X = E * X + D * (Pdyn(:, i) + Pleak);
      end

      %
      % Do not forget the last temperature vector.
      %
      T(:, end) = BT * X + Tamb;
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, hsConfig, hsLine);
  end
end
