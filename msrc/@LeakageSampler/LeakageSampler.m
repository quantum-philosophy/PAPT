classdef LeakageSampler < handle
  properties (Constant)
    Lnom = Leakage.Lnom;
    Ldev = 0.20 * Leakage.Lnom;
  end

  properties (SetAccess = 'private')
    leakage
    pc
    coeff
  end

  methods
    function ls = LeakageSampler(leakage, pc)
      ls.leakage = leakage;
      ls.pc = pc;
    end

    function setup(ls, coeff)
      ls.coeff = coeff;
    end

    function P = perform(ls, samples)
      P = ls.leakage.calculate(ls.calculateTemperature(samples), ...
        ls.Lnom + ls.Ldev .* samples);
    end
  end

  methods (Access = 'private')
    function T = calculateTemperature(ls, rvs)
      T = ls.coeff(:, 1);

      return;

      %
      % Without inner expansions for now.
      %

      if size(ls.coeff, 2) == 1
        T = ls.coeff;
      else
        T = transpose(ls.pc.evaluate(ls.coeff, rvs));
      end
    end
  end
end
