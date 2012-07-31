function x = constructLinearSpace(varargin)
  [ left, right ] = Utils.detectBounds(varargin{:});
  if left < 0 || right > 10^3
    error('Cannot construct a liner space with left = %e and right = %e.', left, right);
  end
  points = max((right - left) / 0.1, 100);
  x = linspace(left, right, points);
end
