classdef PowerChaos < HotSpot.Chaos
  methods
    function this = PowerChaos(varargin)
      this = this@HotSpot.Chaos(varargin{:});
    end
  end

  methods (Access = 'protected')
    function P = solve(this, Pdyn, leakage, rvs)
      [ ~, P ] = solve@HotSpot.Chaos(this, Pdyn, leakage, rvs);
    end
  end
end
