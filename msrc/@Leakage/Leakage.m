classdef Leakage < handle
  properties (Constant)
    %
    % The nominal value of the channel length and its deviation.
    %
    Lnom = 45e-9;
    Ldev = 0.05 * LeakageSampler.Lnom;

    %
    % The following constants are used to construct an instance
    % of the leakage model that produces `PleakPdyn' portion
    % of the given dynamic power at temperature `Tref'.
    %
    PleakPdyn = 2/3;
    Tref = Utils.toKelvin(120);
  end

  properties (SetAccess = 'private')
    %
    % A multivariate polynomial p(L, T) to compute the leakage current.
    %
    compute = Spice.fit('inverter_45nm', [ 2, 2 ]);

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
    function ls = Leakage(Tamb, Pdyn, trace, pc)
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
      P0 = ls.compute(ls.Lnom, ls.Tref);
      ls.alpha = ls.PleakPdyn * mean(Pdyn, 2) ./ P0;

      %
      % Adjust to the dimension of the quadrature.
      %
      ls.alpha = irep(ls.alpha, 1, ls.count);
      ls.Tamb = irep(ls.Tamb, 1, ls.count);
    end

    function advance(ls)
      ls.position = ls.position + 1;
    end

    function P = performAtAmbient(ls, rvs)
      P = ls.alpha .* ls.compute(ls.Lnom + ls.Ldev .* rvs, ls.Tamb);
    end

    function P = performAtCurrent(ls, rvs)
      %
      % With...
      %
      T = ls.pc.evaluate(ls.trace.coeff(:, :, ls.position), rvs);
      %
      % ... or without inner expansions.
      %
      % T = ls.trace.coeff(:, 1, ls.position);
      % T = irep(T, 1, ls.count);
      %
      P = ls.alpha .* ls.compute(ls.Lnom + ls.Ldev .* rvs, T);
    end

    function P = performAtGiven(ls, T, rvs)
      P = ls.alpha .* ls.compute(ls.Lnom + ls.Ldev .* rvs, T);
    end
  end
end
