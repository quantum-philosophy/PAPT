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
    % C^(-1/2)
    Cm12

    % Input and output mapping matrices
    Min
    Mout

    % Min_ = C^(-1/2) * Min
    Min_

    % Mout_ = C^(-1/2) * Mout
    Mout_

    % G_ = C^(-1/2) * (-G) * C^(-1/2)
    G_

    % Eigenvalue decomposition of Gt
    % Gt = V * L * V^T
    L
    V
    VT

    % A = exp(G_ * t) = V * diag(exp(li * dt)) * V^T
    A

    % B = D^(-1) * (exp(G_ * t) - I) * Mint
    %   = V * diag((exp(li * t) - 1) / li) * V^T * Mint
    B
  end

  methods
    function hs = HotSpot(floorplan, hsConfig, hsLine)
      [ C, G, hs.nodes, hs.cores, hs.dt, hs.Tamb ] = ...
        HotSpot.getCoefficients(floorplan, hsConfig, hsLine);

      hs.nodes = size(G, 1);

      hs.Cm12 = diag(sqrt(1 ./ C));

      hs.Min = [ diag(ones(1, hs.cores)); zeros(hs.nodes - hs.cores, hs.cores) ];
      hs.Mout = hs.Min;

      hs.Min_ = hs.Cm12 * hs.Min;
      hs.Mout_ = hs.Cm12 * hs.Mout;

      hs.G_ = Utils.symmetrize(hs.Cm12 * (-G) * hs.Cm12);

      [ V, L ] = eig(hs.G_);

      hs.L = diag(L);
      hs.V = V;
      hs.VT = V';

      hs.A = hs.V * diag(exp(hs.dt * hs.L)) * hs.VT;
      hs.B = hs.V * diag((exp(hs.dt * hs.L) - 1) ./ hs.L) * hs.VT * hs.Min_;
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

      A = hs.A;
      B = hs.B;

      T(:, 1) = A * T0 + B * P(:, 1);

      for i = 2:steps
        T(:, i) = A * T(:, i - 1) + B * P(:, i);
      end

      T = hs.Mout_' * T + hs.Tamb;
    end
  end

  methods (Static)
    [ C, G, nodes, cores, dt, Tamb ] = getCoefficients(floorplan, hsConfig, hsLine);
  end
end
