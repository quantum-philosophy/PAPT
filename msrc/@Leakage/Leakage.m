classdef Leakage < handle
  %
  % Description:
  %
  %   The leakage model is based on
  %
  %   W. Liao, L. He, and K. M. Lepak, "Temperature and Supply Voltage
  %   Aware Performance and Power Modeling at Microarchitecture Level",
  %   IEEE Trans. on Computer-Aided Design of Integrated Circuits and
  %   Systems, July 2005.
  %
  properties (Constant)
    %
    % Scaling coefficients of the average leakage current for 65nm
    %
    A = 1.1432e-12;
    B = 1.0126e-14;
    alpha = 466.4029;
    beta = -1224.74083;
    gamma = 6.28153;
    delta = 6.9094;
    Is = Leakage.calculateMeanIs();
  end

  methods (Static)
    function P = calculate(Ngate, Vdd, T)
      %
      % The average leakage current per gate:
      %
      % Iavg(T, Vdd) = Is(T0, V0) * favg(T, Vdd)
      %
      Iavg = Leakage.Is * Leakage.calculateScaling(T, Vdd);

      %
      % The power leakage for all gates:
      %
      % P = Ngate * Iavg * Vdd
      %
      P = Ngate .* Iavg .* Vdd;
    end
  end

  methods (Static, Access = 'private')
    function f = calculateScaling(T, Vdd)
      %
      % The scaling factor of the leakage current:
      %
      % f(T, Vdd) = A * T^2 * e^((alpha * Vdd + beta)/T) +
      %             B * e^(gamma * Vdd + delta)
      %
      f = Leakage.A .* T.^2 .* exp((Leakage.alpha .* Vdd + Leakage.beta) ./ T) + ...
        Leakage.B .* exp(Leakage.gamma .* Vdd + Leakage.delta);

      % f = 1: Vdd ~= 4.03, the authors say it is not true
    end

    function Is = calculateMeanIs()
      %
      % Some values found in the paper.
      %
      T = Utils.toKelvin([ 100, 100, 80, 80, 60, 60 ]);
      V = [ 0.95, 1.05, 0.95, 1.05, 0.95, 1.05 ];
      Iavg = [ 23.44, 29.56, 19.44, 25.14, 16.00, 21.33 ] * 1e-6;

      Is = mean(Iavg ./ Leakage.calculateScaling(T, V));
    end
  end
end
