classdef Exponent < Leakage.Base
  properties (Constant)
    %
    % The coefficients are found in
    %
    % W. Liao, L. He, and K. M. Lepak, "Temperature and Supply Voltage
    % Aware Performance and Power Modeling at Microarchitecture Level",
    % IEEE Trans. on Computer-Aided Design of Integrated Circuits and
    % Systems, July 2005.
    %
    beta = -(466.4029 * 1 - 1224.74083) / Leakage.Exponent.Lnom;
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
    % Holds the current position in `trace' and is used
    % for sequential sampling.
    %
    position

    %
    % The number of samples that we are going to require at once.
    %
    points
  end

  methods
    function lk = Exponent(Tamb, Pdyn, pca, trace, pc)
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
        lk.points = pc.qd.points;
      end

      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %
      lk.alpha = 1;
      P0 = lk.Tref.^2 .* exp(- lk.beta .* lk.Lnom ./ lk.Tref);
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
      T = lk.Tamb;
      P = lk.alpha .* T.^2 .* exp(- lk.beta .* ...
        (lk.Lnom + lk.pca * rvs) ./ T);
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
      P = lk.alpha .* T.^2 .* exp(- lk.beta .* ...
        (lk.Lnom + lk.pca * rvs) ./ T);
    end

    function P = performAtGiven(lk, T, rvs)
      P = lk.alpha .* T.^2 .* exp(- lk.beta .* ...
        (lk.Lnom + lk.pca * rvs) ./ T);
    end
  end
end
