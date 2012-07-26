function result = findCache(pattern)
  files = dir(Utils.resolvePath('*.mat', 'cache'));

  count = length(files);

  names = {};
  scores = [];

  for i = 1:count
    name = files(i).name;
    M = regexp(name, pattern, 'tokens');
    if ~isempty(M)
      names{end + 1} = name;
      scores(end + 1) = str2num(M{1}{1});
    end
  end

  if isempty(names)
    result = {};
  else
    [ ~, i ] = max(scores);
    result = Utils.resolvePath(names{i}, 'cache');
  end
end
