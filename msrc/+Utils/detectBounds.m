function [ left, right ] = detectBounds(varargin)
  count = length(varargin);

  range = '3sigma';

  if isa(varargin{end}, 'struct')
    count = count - 1;
    options = varargin{end};
    if isfield(options, 'range')
      range = options.range;
    end
  end

  if isa(range, 'char')
    if strcmpi(range, 'unbounded')
      range = @(~, ~) [ -Inf, Inf ];
    else
      tokens = regexp(range, '^(.+)sigma$', 'tokens');
      if isempty(tokens), error('The range in unknown.'); end
      times = str2num(tokens{1}{1});
      range = @(mu, sigma) [ mu - times * sigma, mu + times * sigma ];
    end
  end

  left = zeros(count, 1);
  right = zeros(count, 1);

  for i = 1:count
    one = varargin{i};

    mn = min(one);
    mx = max(one);

    mu = mean(one);
    sigma = sqrt(var(one));

    r = range(mu, sigma);

    left(i) = max(mn, r(1));
    right(i) = min(mx, r(2));
  end

  left = min(left);
  right = max(right);
end
