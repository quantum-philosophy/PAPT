function localError = compare(mcRaw, pcRaw, varargin)
  if length(varargin) == 1 && isa(varargin{1}, 'struct')
    options = varargin{1};
  else
    options = struct(varargin{:});
  end

  [ ~, ddim ] = size(mcRaw);

  draw = Utils.extract(options, 'draw', false);

  if draw, figure; end

  localError = zeros(1, ddim);

  for i = 1:ddim
    if draw, p = subplot(1, ddim, i); end

    mcraw = mcRaw(:, i);
    pcraw = pcRaw(:, i);

    x = Utils.constructLinearSpace(mcraw, pcraw, options);

    [ mcx, mcData ] = process(x, mcraw, options);
    [ pcx, pcData ] = process(x, pcraw, options);

    assert(nnz(mcx - pcx) == 0, 'The X vectors are invalid.');

    localError(i) = Utils.NRMSE(mcData, pcData);

    if draw
      paint(mcx, mcData, pcData, options);
      title(sprintf('NRMSE %.4e', localError(i)));
      labels = Utils.extract(options, 'labels', { 'MC', 'PC' });
      legend(labels{:});
    end
  end
end

function [ x, data ] = process(x, raw, options)
  method = Utils.extract(options, 'method', 'smooth');

  switch method
  case 'smooth'
    [ x, data ] = processKernelDensity(x, raw, options);
  case 'histogram'
    [ x, data ] = processHistogram(x, raw, options);
  case 'piecewise'
    [ x, data ] = processHistogram(x, raw, options);
  otherwise
    error('The method is unknown.');
  end
end

function paint(x, mcData, pcData, options)
  method = Utils.extract(options, 'method', 'smooth');

  switch method
  case 'smooth'
    paintLines(x, mcData, pcData, options);
  case 'histogram'
    paintHistogram(x, mcData, pcData, options);
  case 'piecewise'
    paintLines(x, mcData, pcData, options);
  otherwise
    error('The method is unknown.');
  end
end

function [ x, data ] = processHistogram(x, raw, options)
  x = x(:);
  raw = raw(:);

  dx = x(2:end) - x(1:end - 1);
  dx = [ dx(1); dx; dx(end) ];

  data = histc(raw, [ -Inf; x; Inf ]);
  data = data(1:end - 1);

  data = data ./ dx / sum(data);
  x = [ x(1) - dx(1); x ];

  switch Utils.extract(options, 'function', 'pdf');
  case 'pdf'
  case 'cdf'
    data = cumsum(data);
    data = data / data(end);
  otherwise
    error('The function is unknown.');
  end
end

function [ x, data ] = processKernelDensity(x, raw, options)
  data = ksdensity(raw, x, ...
    'function', Utils.extract(options, 'function', 'pdf'));
end

function paintHistogram(x, mcData, pcData, options)
  hold on;
  drawBar(x, mcData, Utils.pickColor(1));
  drawBar(x, pcData, Utils.pickColor(2));
end

function drawBar(x, data, color)
  hbar = bar(x, data, 'FaceColor', color, 'Edgecolor', color);
  hpatch = findobj(gca, 'Type', 'patch');
  set(hpatch, 'FaceAlpha', 0.75);
end

function paintLines(x, mcData, pcData, options)
  line(x, mcData, 'Color', Utils.pickColor(1));
  line(x, pcData, 'Color', Utils.pickColor(2));
end
