classdef LeakageSampler < handle
  properties (Constant)
    %
    % The nominal value of the channel length and its deviation.
    %
    Lnom = 45e-9;
    Ldev = 0.05 * LeakageSampler.Lnom;

    %
    % The coefficients are found in
    %
    % W. Liao, L. He, and K. M. Lepak, "Temperature and Supply Voltage
    % Aware Performance and Power Modeling at Microarchitecture Level",
    % IEEE Trans. on Computer-Aided Design of Integrated Circuits and
    % Systems, July 2005.
    %
    beta = -(466.4029 * 1 - 1224.74083) / LeakageSampler.Lnom;

    %
    % The following constants are used to construct an instance
    % of the leakage model that produces `PleakPdyn' portion
    % of the given dynamic power at temperature `Tref'.
    %
    PleakPdyn = 1.0;
    Tref = Utils.toKelvin(80);
  end

  properties (SetAccess = 'private')
    %
    % The temperature of the ambience; it is used for the very
    % first sampling round.
    %
    Tamb

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
    % One of the leakage parameters.
    %
    alpha

    %
    % The number of cores that we are going to feed with leakage.
    %
    cores

    %
    % The number of samples that we are going to require at once.
    %
    count
  end

  methods
    function ls = LeakageSampler(Tamb, Pdyn, trace, pc)
      if nargin < 3
        %
        % We are working with Monte-Carlo, one sample at a time.
        %
        ls.count = 1;
      else
        ls.pc = pc;
        ls.trace = trace;
        ls.position = 0;
        ls.count = pc.cq.count;
      end

      ls.cores = size(Pdyn, 1);

      ls.Tamb = ones(ls.cores, 1) * Tamb;

      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %
      ls.alpha = 1;
      P0 = ls.Tref.^2 .* exp(- ls.beta .* ls.Lnom ./ ls.Tref);
      ls.alpha = ls.PleakPdyn * mean(Pdyn, 2) ./ P0;

      %
      % Adjust to the dimension of the quadrature.
      %
      ls.alpha = repmat(ls.alpha, 1, ls.count);
      ls.Tamb = repmat(ls.Tamb, 1, ls.count);
    end

    function advance(ls)
      ls.position = ls.position + 1;
    end

    function P = performAtAmbient(ls, rvs)
      T = ls.Tamb;
      P = ls.alpha .* T.^2 .* exp(- ls.beta .* (ls.Lnom + ls.Ldev .* rvs) ./ T);
    end

    function P = performAtCurrent(ls, rvs)
      %
      % With or without inner expansions.
      %
      T = ls.trace.coeff(:, 1, ls.position);
      % T = ls.pc.evaluate(ls.trace.coeff(:, :, ls.position), rvs);

      T = repmat(T, 1, ls.count);
      P = ls.alpha .* T.^2 .* exp(- ls.beta .* (ls.Lnom + ls.Ldev .* rvs) ./ T);
    end

    function P = performAtGiven(ls, T, rvs)
      P = ls.alpha .* T.^2 .* exp(- ls.beta .* (ls.Lnom + ls.Ldev .* rvs) ./ T);
    end
  end
end
