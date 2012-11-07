function C = construct(this, floorplan, options)
  D = dlmread(floorplan, '', 0, 1);

  W = D(:, 1);
  H = D(:, 2);
  X = D(:, 3);
  Y = D(:, 4);

  dieW = max(X + W);
  dieH = max(Y + H);

  processorX = X + W / 2 - dieW / 2;
  processorY = Y + H / 2 - dieH / 2;

  domain = max(dieW, dieH) / 2;

  kl = KarhunenLoeve.Exponential( ...
    'domainBoundary', domain, ...
    'correlationLength', domain, ...
    'threshold', 1 - options.get('threshold', this.threshold));

  processorCount = size(D, 1);
  dimensionCount = kl.dimension;

  C = zeros(processorCount, dimensionCount);

  for i = 1:processorCount
    distance = sqrt(processorX(i)^2 + processorY(i)^2);
    for j = 1:dimensionCount
      C(i, j) = sqrt(kl.values(j)) * kl.functions{j}(distance);
    end
  end
end
