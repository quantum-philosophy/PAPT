classdef Discrete < ProcessVariation.Base
  methods
    function this = Discrete(varargin)
      this = this@ProcessVariation.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    mapping = construct(this, floorplan, options)
  end
end
