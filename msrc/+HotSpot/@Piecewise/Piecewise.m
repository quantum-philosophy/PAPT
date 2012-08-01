classdef Piecewise < HotSpot.Analytic
  properties (Access = 'protected')
    pieces

    Trange
    Icoeff

    B
    dCm1
    dCm12

    %
    % For the Woodbury formula.
    %
    Am1
    Am1U
    VAm1
    VAm1U
  end

  methods
    function hs = Piecewise(varargin)
      hs = hs@HotSpot.Base(varargin{:});

      [ hs.Trange, hs.Icoeff ] = Spice.fitPiecewiseLinear(...
        Leakage.Polynomial.compute, Leakage.Polynomial.Lnom);

      hs.pieces = size(hs.Trange, 1);

      hs.B = hs.BT';
      hs.dCm1 = 1 ./ hs.C;
      hs.dCm12 = sqrt(1 ./ hs.C);

      %
      % Prepare the Woodbury formula.
      %
      U = [ diag(ones(hs.cores, 1)); zeros(hs.nodes - hs.cores, hs.cores) ];
      V = U';

      hs.Am1 = hs.Q * diag(1 ./ hs.L) * hs.QT;
      hs.Am1U = hs.Am1 * U;
      hs.VAm1 = V * hs.Am1;
      hs.VAm1U = hs.VAm1 * U;
    end

    function [ T, Pleak ] = solve(hs, Pdyn, rvs)
      [ cores, steps ] = size(Pdyn);
      assert(cores == hs.cores, 'The power profile is invalid.')

      %
      % The goals.
      %
      T = zeros(cores, steps);
      Pleak = zeros(size(Pdyn));

      %
      % General shortcuts.
      %
      nodes = hs.nodes;
      E = hs.E;
      D = hs.D;
      BT = hs.BT;
      Tamb = hs.Tamb;

      dt = hs.dt;
      B = hs.B;

      pieces = hs.pieces;

      dCm1 = hs.dCm1;
      dCm12 = hs.dCm12;

      Trange = hs.Trange;
      Icoeff = hs.Icoeff;

      Am1 = hs.Am1;
      Am1U = hs.Am1U;
      VAm1 = hs.VAm1;
      VAm1U = hs.VAm1U;

      I = diag(ones(nodes, 1));

      %
      % Initialize the leakage model.
      %
      leak = Leakage.Polynomial(Tamb, Pdyn, hs.pca);
      alpha = leak.alpha;
      Kcoeff = zeros(pieces, cores);
      Mcoeff = zeros(pieces, cores);
      Lcoeff = zeros(pieces, cores);
      for i = 1:cores
        Kcoeff(:, i) = alpha(i) *            Icoeff(:, 2);
        Mcoeff(:, i) = alpha(i) *            Icoeff(:, 1);
        Lcoeff(:, i) = alpha(i) * dCm1(i)  * Icoeff(:, 2);
      end

      %
      % Determine the scale coefficients.
      %
      k = min(find(Trange(:, 2) > Tamb));
      if isempty(k), k = pieces; end;

      K = zeros(cores, 1);
      M = zeros(cores, 1);
      L = zeros(nodes, 1);

      K(1:cores) = Kcoeff(k, :);
      M(1:cores) = Mcoeff(k, :);
      L(1:cores) = Lcoeff(k, :);

      %
      % Update the coefficient matrices.
      %
      % E_ = E * diag(exp(L * dt));
      [ V, C ] = eig(hs.A + diag(L));
      C = diag(C);
      E_ = V * diag(exp(C * dt)) * V';
      hs.E = hs.Q * diag(exp(hs.dt * hs.L)) * hs.QT;

      % D_ = (Am1 - Am1U * inv(diag(1 ./ L(1:cores)) + VAm1U) * VAm1) * (E_ - I) * B;
      D_ = inv(hs.A + diag(L)) * (E_ - I) * B;

      %
      % Compute the leakage power.
      %
      % Pleak(:, 1) = leak.performAtAmbient(rvs);
      Pleak(:, 1) = M;

      %
      % Compute the state vector.
      %
      X = D_ * (Pdyn(:, 1) + K * Tamb + M);

      for i = 2:steps
        T(:, i - 1) = BT * X + Tamb;

        %
        % Determine the scale coefficients.
        %
        for j = 1:cores
          k = min(find(Trange(:, 2) > T(j, i - 1)));
          if isempty(k), k = pieces; end;
          K(j) = Kcoeff(k, j);
          M(j) = Mcoeff(k, j);
          L(j) = Lcoeff(k, j);
        end

        %
        % Update the coefficient matrices.
        %
        % E_ = E * diag(exp(L * dt));
        [ V, C ] = eig(hs.A + diag(L));
        C = diag(C);
        E_ = V * diag(exp(C * dt)) * V';
        % D_ = (Am1 - Am1U * inv(diag(1 ./ L(1:cores)) + VAm1U) * VAm1) * (E_ - I) * B;
        D_ = inv(hs.A + diag(L)) * (E_ - I) * B;

        %
        % Compute the leakage power.
        %
        % Pleak(:, i) = leak.performAtGiven(rvs, T(:, i - 1));
        Pleak(:, i) = M;

        %
        % Compute the state vector.
        %
        X = E_ * X + D_ * (Pdyn(:, i) + K * Tamb + M);
      end

      %
      % Do not forget the last temperature vector.
      %
      T(:, end) = BT * X + Tamb;
    end
  end
end
