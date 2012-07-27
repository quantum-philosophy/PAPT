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

  methods
    function lk = Exponent(Tamb, Pdyn, pca, points)
      [ cores, steps ] = size(Pdyn);

      lk = lk@Leakage.Base(Tamb, cores, pca);

      if nargin < 4
        %
        % One sample at a time (Monte-Carlo).
        %
        points = 1;
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
      lk.alpha = irep(lk.alpha, 1, points);
      lk.Tamb = irep(lk.Tamb, 1, points);
    end

    function P = performAtAmbient(lk, rvs)
      T = lk.Tamb;
      P = lk.alpha .* T.^2 .* exp(- lk.beta .* ...
        (lk.Lnom + lk.pca * rvs) ./ T);
    end

    function P = performAtGiven(lk, rvs, T)
      P = lk.alpha .* T.^2 .* exp(- lk.beta .* ...
        (lk.Lnom + lk.pca * rvs) ./ T);
    end
  end
end
