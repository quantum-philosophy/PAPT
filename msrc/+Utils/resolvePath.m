function fullpath = path(file, type)
  if file(1) == '#'
    fullpath = [ pwd, '/', file(2:end) ];
  else
    if nargin < 2, type = 'working'; end
    fullpath = [ Constants.([ type, 'Directory' ]), '/', file ];
  end
end
