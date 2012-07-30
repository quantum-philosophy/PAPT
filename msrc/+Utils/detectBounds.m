function [ left, right ] = detectBounds(varargin)
  count = length(varargin);

  left = zeros(count, 1);
  right = zeros(count, 1);

  for i = 1:count
    one = varargin{i};

    mn = min(one);
    mx = max(one);

    exp = mean(one);
    std = sqrt(var(one));

    left(i) = max(mn, exp - 3 * std);
    right(i) = min(mx, exp + 3 * std);
  end

  left = min(left);
  right = max(right);
end
