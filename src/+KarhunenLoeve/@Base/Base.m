classdef Base < handle
  properties (SetAccess = 'private')
    dimension
    domainBoundary
    correlationLength

    values
    functions
  end

  methods
    function this = Base(varargin)
      options = Options(varargin{:});
      this.initialize(options);
    end
  end

  methods (Abstract, Access = 'protected')
    [ values, functions ] = construct(this, options)
  end

  methods (Access = 'private')
    function initialize(this, options)
      this.domainBoundary = options.domainBoundary;
      this.correlationLength = options.correlationLength;
      [ this.values, this.functions ] = this.construct(options);
      this.dimension = length(this.values);
    end
  end
end
