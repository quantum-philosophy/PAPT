function x = constructLinearSpace(varargin)
  [ left, right ] = Stats.detectBounds(varargin{:});

  if isa(varargin{end}, 'struct')
    options = varargin{end};
  else
    options = struct();
  end

  if isfield(options, 'points')
    points = options.points;
  else
    points = max((right - left) / 0.1, 100);
  end

  x = linspace(left, right, points);
end
