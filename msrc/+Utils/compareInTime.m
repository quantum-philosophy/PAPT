function [ globalError, localError ] = compareInTime(mcRaw, pcRaw, varargin)
  if length(varargin) == 1 && isa(varargin{1}, 'struct')
    options = varargin{1};
  else
    options = struct(varargin{:});
  end

  [ ~, ddim, tdim ] = size(mcRaw);

  draw = Utils.extract(options, 'draw', false);
  options.draw = false;

  h = ibar('Comparison of PDFs: step %d out of %d.', tdim);

  localError = zeros(ddim, tdim);

  for i = 1:tdim
    localError(:, i) = Utils.compare(mcRaw(:, :, i), pcRaw(:, :, i), options);
    increase(h);
  end

  globalError = sqrt(sum(localError(:) .^ 2) / (ddim * tdim));

  if ~draw, return; end

  figure;

  labels = {};
  time = 0:(tdim - 1);

  for i = 1:ddim
    labels{end + 1} = sprintf('Dimension %d', i);
    line(time, localError(i, :) * 100, 'Color', Utils.pickColor(i));
  end

  title('Evolution of Error');
  ylabel('NRMSE, %');
  legend(labels{:});
end
