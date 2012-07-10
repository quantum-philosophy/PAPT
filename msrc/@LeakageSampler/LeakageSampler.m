classdef LeakageSampler < Leakage
  properties (SetAccess = 'private')
    cores

    Vdd
    Ngate
  end

  methods
    function ls = LeakageSampler(P, T)
      %
      % Fit the leakage coefficients to produce the leakage power `P'
      % at the temperature level `T'.
      %

      ls.cores = length(P);

      if nargin < 2, T = 100; end
      if length(T) == 1, T = ones(ls.cores, 1) * T; end

      %
      % Let us fix Vdd and find Ngate.
      %
      ls.Vdd = ones(ls.cores, 1); % V
      ls.Ngate = ones(ls.cores, 1);

      P0 = Leakage.calculate(ls.Ngate, ls.Vdd, T);

      ls.Ngate = floor(P ./ P0);
    end

    function setup(ps, varargin)
    end

    function P = sample(ls, values)
      P = values;
    end
  end
end
