function [ data, options ] = extract(varargin)
  count = length(varargin);

  i = 1;
  while i <= count && isa(varargin{i}, 'double'); i = i + 1; end
  data = varargin(1:(i - 1));

  options = Options(varargin{i:end});
end
