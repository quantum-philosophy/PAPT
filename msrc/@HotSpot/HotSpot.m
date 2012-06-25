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
    % C^{-1/2}
    Cm12

    % Eigenvalue decomposition
    % C^(-1/2) * (-G) * C^(-1/2) = V * L * V^T
    L
    V
    VT

    % A = V * diag(exp(li * dt)) * V^T
    A

    % B = D^(-1) (exp(D * t) - I) C^(-1/2)
    %   = V * diag((exp(li * t) - 1) / li) * V^T * C^(-1/2)
    B
  end

  methods
    function hs = HotSpot(floorplan, hsConfig, hsLine)
      [ C, G, hs.nodes, hs.cores, hs.dt, hs.Tamb ] = ...
        HotSpot.getCoefficients(floorplan, hsConfig, hsLine);

      hs.nodes = size(G, 1);

      hs.Cm12 = diag(sqrt(1 ./ C));

      [ V, L ] = eig(Utils.symmetrize(hs.Cm12 * (-G) * hs.Cm12));

      hs.L = diag(L);
      hs.V = V;
      hs.VT = V';

      hs.A = hs.V * diag(exp(hs.dt * hs.L)) * hs.VT;
      hs.B = hs.V * diag((exp(hs.dt * hs.L) - 1) ./ hs.L) * hs.VT * hs.Cm12;
    end

    function T = solve(hs, P, T0)
      [ cores, steps ] = size(P);
      nodes = hs.nodes;
      Tamb = hs.Tamb;

      A = hs.A;
      B = hs.B;
      Cm12 = hs.Cm12;

      if cores ~= hs.cores, error('The power profile is invalid.'); end

      if nargin < 3
        T0 = zeros(nodes, 1);
      else
        T0 = Cm12.^(-1) * (T0 - Tamb);
      end

      T = zeros(nodes, steps);

      padding = zeros(nodes - cores, 1);

      T(:, 1) = A * T0 + B * [ P(:, 1); padding ];

      for i = 2:steps
        T(:, i) = A * T(:, i - 1) + B * [ P(:, i); padding ];
        T(:, i - 1) = Cm12 * T(:, i - 1) + Tamb;
      end

      T(:, steps) = Cm12 * T(:, steps) + Tamb;

      T = T(1:cores, :);
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, hsConfig, hsLine);
  end
end
