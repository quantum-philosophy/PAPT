classdef Polynomial < Leakage.Base
  properties (Constant)
    %
    % A multivariate polynomial p(L, T) to compute the leakage current.
    %
    compute = Spice.fitExponentPolynomial( ...
      'inverter_45nm', [ 1, 2 ], [ 1, 0.7, 0; 1, 1, 1 ]);
  end

  methods
    function lk = Polynomial(Tamb, Pdyn, pca, points)
      [ cores, steps ] = size(Pdyn);

      lk = lk@Leakage.Base(Tamb, cores, pca);

      if nargin < 4
        %
        % One sample at a time (Monte Carlo).
        %
        points = 1;
      end

      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %
      P0 = lk.compute(lk.Lnom, lk.Tref);
      lk.alpha = lk.PleakPdyn * mean(mean(Pdyn)) ./ P0;

      %
      % Adjust to the dimension of the quadrature.
      %
      lk.Tamb = irep(lk.Tamb, 1, points);
    end

    function P = performAtAmbient(lk, rvs)
      P = lk.alpha .* lk.compute(lk.Lnom + lk.pca * rvs, lk.Tamb);
    end

    function P = performAtGiven(lk, rvs, T)
      P = lk.alpha .* lk.compute(lk.Lnom + lk.pca * rvs, T);
    end
  end
end
