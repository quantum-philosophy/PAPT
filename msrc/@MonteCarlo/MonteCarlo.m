classdef MonteCarlo < handle
  properties (Constant)
    Lnom = Leakage.Lnom;
    Ldev = 0.20 * Leakage.Lnom;
  end

  properties (Access = 'private')
    hotspot
    Pdyn
    leakage
  end

  methods
    function mc = MonteCarlo(hotspot, Pdyn, leakage)
      mc.hotspot = hotspot;
      mc.Pdyn = Pdyn;
      mc.leakage = leakage;
    end

    function T = evaluate(mc, rvs)
      mc.leakage.adjust(mc.Lnom + mc.Ldev .* rvs);
      T = mc.hotspot.solve(mc.Pdyn, mc.leakage);
    end
  end

  methods (Static)
    [ E, C, out ] = perform(f, dims, samples);
    [ E, C ] = perform3D(f, dims, samples);
  end
end
