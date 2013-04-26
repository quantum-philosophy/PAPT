classdef PowerChaos < HotSpot.Chaos
  methods
    function this = PowerChaos(varargin)
      this = this@HotSpot.Chaos(varargin{:});
    end
  end

  methods (Access = 'protected')
    function P = solve(this, varargin);
      [ ~, P ] = solve@HotSpot.Chaos(this, varargin{:});
    end
  end
end
