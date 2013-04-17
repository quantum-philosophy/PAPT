classdef Exponential < KarhunenLoeve.Base
  methods
    function this = Exponential(varargin)
      this = this@KarhunenLoeve.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    [ values, functions ] = construct(this, options)
  end
end
