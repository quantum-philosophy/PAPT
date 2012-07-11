function f = toFunction(p, varargin)
  %
  % Description:
  %
  %   Converts the given polynomial `p' to a function with multiple
  %   vector-valued inputs formed according to the given sets of variables.
  %
  % Example:
  %
  % >> toFunction(a0*x2^2 + a0 + a1*x1 + x3^3, [ x1 x2 x3 ], [ a0 a1 ])
  % ans =
  % @(y1, y2) y2(1) .* y1(2).^2 + y2(1) + y2(2) .* y1(1) + y1(3).^3
  %

  count = length(varargin);

  y = {};
  t = {};

  pending = false;

  for i = 1:count
    if ~ischar(varargin{i})
      y{end + 1} = varargin{i};
      if pending
        t{end + 1} = 'y%d(%d)';
      end
      pending = true;
    else
      orientation = lower(varargin{i});
      switch (orientation)
      case 'columns'
        t{end + 1} = 'y%d(:,%d)';
      case 'rows'
        t{end + 1} = 'y%d(%d,:)';
      otherwise
        error('The specified orientation is unknown.');
      end
      pending = false;
    end
  end

  if pending
    t{end + 1} = 'y%d(%d)';
  end

  count = length(y);
  args = 'y1';

  for i = 2:count
    args = [ args, sprintf(',y%d', i) ];
  end

  if isa(p, 'sympoly')
    s = char(p);
    s = regexprep(s, '\^', '.^');
    s = regexprep(s, '\*', '.*');
    s = regexprep(s, '\/', './');
  else
    f = matlabFunction(p);
    s = func2str(f);
    s = regexprep(s, '@\([^)]+\)', '');
  end

  for i = 1:count
    y0 = y{i};
    for j = 1:length(y0);
      s = regexprep(s, [ '\<', char(y0(j)), '\>' ], sprintf(t{i}, i, j));
    end
  end

  f = str2func([ '@(', args, ')', s ]);
end
