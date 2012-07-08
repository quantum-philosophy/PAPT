classdef HotSpot < handle
  properties (SetAccess = public)
    % Number of thermal nodes
    nodes

    % Number of cores
    cores

    % Sampling interval
    dt

    % Ambient temperature
    Tamb
  end

  properties (SetAccess = private)
    % Original system:
    % C * dZ/dt + G * Z = M * P
    % T = M^T * Z + T_amb

    % Transformed system:
    % dX/dt = A * X + B * P
    % T = B^T * X + T_amb

    % C^(-1/2)
    Cm12

    % A = C^(-1/2) * (-G) * C^(-1/2)
    A

    % B^T = (C^(-1/2) * M)^T
    BT

    % Eigenvalue decomposition of Gt
    % A = Q * L * Q^T
    L
    Q
    QT

    % E = exp(A * t) = Q * diag(exp(li * dt)) * Q^T
    E

    % D = A^(-1) * (exp(A * t) - I) * B
    %   = Q * diag((exp(li * t) - 1) / li) * Q^T * B
    D
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
    end

    function T = solve(hs, P, T0)
      [ cores, steps ] = size(P);
      nodes = hs.nodes;

      if cores ~= hs.cores, error('The power profile is invalid.'); end

      if nargin < 3
        T0 = zeros(nodes, 1);
      else
        T0 = hs.Cm12.^(-1) * (T0 - Tamb);
      end

      T = zeros(nodes, steps);

      E = hs.E;
      D = hs.D;

      T(:, 1) = E * T0 + D * P(:, 1);

      for i = 2:steps
        T(:, i) = E * T(:, i - 1) + D * P(:, i);
      end

      T = hs.BT * T + hs.Tamb;
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, hsConfig, hsLine);
  end
end
