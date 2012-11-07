classdef Base < handle
  properties (Constant)
    Lnom = LeakagePower.Lnom;
    Ldev = LeakagePower.Lnom * 0.05;
  end

  properties (SetAccess = 'protected')
    Lmap
    dimension
  end

  methods
    function this = Base(floorplan, varargin)
      options = Options(varargin{:});
      this.Lmap = this.construct(floorplan, options);
      this.dimension = size(this.Lmap, 2);
    end
  end

  methods (Abstract, Access = 'protected')
    mapping = construct(this, floorplan, options)
  end
end
