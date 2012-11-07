classdef Discrete < ProcessVariation.Base
  properties (Constant)
    %
    % The portion of the information that is to be preserved.
    %
    threshold = 0.99;

    %
    % The contribution of the global variations.
    %
    globalPortion = 0.5;
  end

  methods
    function this = Discrete(varargin)
      this = this@ProcessVariation.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function mapping = construct(this, floorplan, options)
      C = computeCorrelation(floorplan);
      P = performPCA(C, options.get('reduction', 'adjustable'), this.threshold);

      portion = this.globalPortion;
      processorCount = size(C, 1);

      mapping = [ sqrt(1 - portion) * P, sqrt(portion) * ones(processorCount, 1) ];
    end
  end
end
