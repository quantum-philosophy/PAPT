classdef Continuous < ProcessVariation.Base
  properties (Constant)
    %
    % The portion of the information that is to be preserved.
    %
    threshold = 0.95
  end

  methods
    function this = Continuous(varargin)
      this = this@ProcessVariation.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function mapping = construct(this, floorplan, options)
      mapping = computeCorrelation(floorplan, this.threshold);
    end
  end
end
