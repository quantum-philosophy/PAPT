function [ x, dx ] = constructLinearSpace(varargin)
  [ left, right ] = Utils.detectBounds(varargin{:});
  points = max((right - left) / 0.1, 100);
  dx = (right - left) / (points - 1);
  x = linspace(left, right, points);
end
