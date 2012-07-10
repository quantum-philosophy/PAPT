classdef Leakage < handle
  properties (Constant)
    %
    % The coefficients are found in
    %
    % W. Liao, L. He, and K. M. Lepak, "Temperature and Supply Voltage
    % Aware Performance and Power Modeling at Microarchitecture Level",
    % IEEE Trans. on Computer-Aided Design of Integrated Circuits and
    % Systems, July 2005.
    %
    Lnom = 45e-9;
    beta = -(466.4029 * 1 -1224.74083) / Leakage.Lnom;

    %
    % The following constants are used to construct an instance
    % of the leakage model that produces `relativeToDynamic' portion
    % of the given dynamic power at temperature `referenceTemperature'.
    %
    relativeToDynamic = 1.0;
    referenceTemperature = Utils.toKelvin(80);
  end

  properties (SetAccess = 'private')
    L
    alpha
  end

  methods
    function l = Leakage(P, T)
      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %

      cores = length(P);

      if nargin < 2, T = Utils.toKelvin(100); end
      if length(T) == 1, T = ones(cores, 1) * T; end

      l.L = ones(cores, 1) * Leakage.Lnom;
      l.alpha = ones(cores, 1);

      P0 = l.calculate(T);

      l.alpha = P ./ P0;
    end

    function P = calculate(l, T, L)
      if nargin < 3, L = l.L; end
      P = l.alpha .* T.^2 .* exp(- l.beta .* L ./ T);
    end
  end

  methods (Static)
    function leakage = constructBasedOnDynamic(Pdyn, R, T)
      %
      % Description:
      %
      %   Construct the leakage power model and adjust it to produce
      %   a curtain about of power relative to the given dynamic power.
      %
      %

      if nargin < 2, R = Leakage.relativeToDynamic; end
      if nargin < 3, T = Leakage.referenceTemperature; end

      leakage = Leakage(R * mean(Pdyn, 2), T);
    end
  end
end
