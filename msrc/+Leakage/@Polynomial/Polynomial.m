classdef Polynomial < Leakage.Base
  properties (Constant)
    %
    % A multivariate polynomial p(L, T) to compute the leakage current.
    %
    % compute = Spice.fitPolynomial('inverter_45nm', [ 2, 2 ]);
    compute = Spice.fitExponentPolynomial('inverter_45nm', [ 2, 2 ], 0.7);
  end

  properties (SetAccess = 'private')
    %
    % The storage of all PC coefficients at all moments of time.
    %
    trace

    %
    % The polynomial chaos used to calculate temperature.
    %
    pc

    %
    % The number of samples that we are going to require at once.
    %
    points
  end

  properties
    %
    % Holds the current position in `trace' and is used
    % for sequential sampling.
    %
    % NOTE: It is public since we want to optimize the core.
    %
    position
  end

  methods
    function lk = Polynomial(Tamb, Pdyn, pca, trace, pc)
      [ cores, steps ] = size(Pdyn);

      lk = lk@Leakage.Base(Tamb, cores, pca);

      if nargin < 4
        %
        % One sample at a time (Monte-Carlo).
        %
        lk.points = 1;
      else
        %
        % A bunch of samples at a time (PC).
        %
        lk.trace = trace;
        lk.pc = pc;
        lk.position = 0;
        lk.points = pc.gq.points;
      end

      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %
      P0 = lk.compute(lk.Lnom, lk.Tref);
      lk.alpha = lk.PleakPdyn * mean(Pdyn, 2) ./ P0;

      %
      % Adjust to the dimension of the quadrature.
      %
      lk.alpha = irep(lk.alpha, 1, lk.points);
      lk.Tamb = irep(lk.Tamb, 1, lk.points);
    end

    function advance(lk)
      lk.position = lk.position + 1;
    end

    function P = performAtAmbient(lk, rvs)
      P = lk.alpha .* lk.compute(lk.Lnom + lk.pca * rvs, lk.Tamb);
    end

    function P = performAtCurrent(lk, rvs)
      %
      % With...
      %
      T = lk.pc.evaluate(lk.trace.coeff(:, :, lk.position), rvs);
      %
      % ... or without inner expansions.
      %
      % T = lk.trace.coeff(:, 1, lk.position);
      % T = irep(T, 1, lk.points);
      %
      P = lk.alpha .* lk.compute(lk.Lnom + lk.pca * rvs, T);
    end

    function P = performAtGiven(lk, T, rvs)
      P = lk.alpha .* lk.compute(lk.Lnom + lk.pca * rvs, T);
    end
  end
end
