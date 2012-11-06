classdef KarhunenLoeve < handle
  properties (SetAccess = 'private')
    dimension
    domainBoundary
    correlationLength

    values
    functions
  end

  methods
    function this = KarhunenLoeve(varargin)
      options = Options(varargin{:});
      this.initialize(options);
    end
  end

  methods (Access = 'private')
    [ values, functions ] = construct(this, options)

    function initialize(this, options)
      this.dimension = options.dimension;
      this.domainBoundary = options.domainBoundary;
      this.correlationLength = options.correlationLength;
      [ this.values, this.functions ] = this.construct(options);
    end
  end
end
