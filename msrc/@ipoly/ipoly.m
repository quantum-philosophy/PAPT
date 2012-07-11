classdef ipoly < sympoly
  methods
    function p = ipoly(varargin)
      p = p@sympoly(varargin{:});
    end

    function np = subs(p, vars, vals)
      n = length(vars);
      np = p;
      for i = 1:n
        np = subs@sympoly(np, vars(i), vals(i));
      end
    end
  end
end
