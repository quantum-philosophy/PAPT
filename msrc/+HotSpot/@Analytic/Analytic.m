classdef Analytic < HotSpot.Base
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
  end

  methods
    function hs = Analytic(floorplan, config, line)
      hs = hs@HotSpot.Base(floorplan, config, line);

      hs.Cm12 = diag(sqrt(1 ./ hs.C));

      M = [ diag(ones(1, hs.cores)); zeros(hs.nodes - hs.cores, hs.cores) ];

      hs.A = Utils.symmetrize(hs.Cm12 * (-hs.G) * hs.Cm12);

      B = hs.Cm12 * M;
      hs.BT = B';

      [ Q, L ] = eig(hs.A);

      hs.L = diag(L);
      hs.Q = Q;
      hs.QT = Q';

      hs.E = hs.Q * diag(exp(hs.dt * hs.L)) * hs.QT;
      hs.D = hs.Q * diag((exp(hs.dt * hs.L) - 1) ./ hs.L) * hs.QT * B;
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
      leak = Leakage.Polynomial(Tamb, Pdyn, hs.pca);

      T = zeros(cores, steps);

      Pleak = leak.performAtAmbient(rvs);
      X = D * (Pdyn(:, 1) + Pleak);

      for i = 2:steps
        T(:, i - 1) = BT * X + Tamb;
        Pleak = leak.performAtGiven(T(:, i - 1), rvs);
        X = E * X + D * (Pdyn(:, i) + Pleak);
      end

      %
      % Do not forget the last temperature vector.
      %
      T(:, end) = BT * X + Tamb;
    end
  end
end
