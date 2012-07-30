function x = constructLinearSpace(varargin)
  [ left, right ] = Utils.detectBounds(varargin{:});
  points = max((right - left) / 0.1, 100);
  x = linspace(left, right, points);
end
