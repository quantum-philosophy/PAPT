function [ left, right ] = detectBounds(varargin)
  count = length(varargin);

  if ischar(varargin{end})
    count = count - 1;
    method = lower(varargin{end});
  else
    method = 'bounded';
  end

  left = zeros(count, 1);
  right = zeros(count, 1);

  for i = 1:count
    one = varargin{i};

    mn = min(one);
    mx = max(one);

    switch method
    case 'bounded'
      exp = mean(one);
      std = sqrt(var(one));

      left(i) = max(mn, exp - 3 * std);
      right(i) = min(mx, exp + 3 * std);
    case 'unbounded'
      left(i) = mn;
      right(i) = mx;
    otherwise
      error('The method is unknown.');
    end
  end

  left = min(left);
  right = max(right);
end
