function mapping = construct(this, options)
  threshold = options.get('threshold', 0.99);
  globalPortion = options.get('globalPortion', 0.5);

  D = dlmread(options.floorplan, '', 0, 1);

  W = D(:, 1);
  H = D(:, 2);
  X = D(:, 3);
  Y = D(:, 4);

  processorCount = size(D, 1);

  dieW = max(X + W);
  dieH = max(Y + H);

  dieX = dieW / 2;
  dieY = dieH / 2;

  processorX = X + W / 2;
  processorY = Y + H / 2;

  correlationLength = max(dieW, dieH) / 2;

  distance = sqrt((dieX - processorX).^2 + (dieY - processorY).^2);

  C = eye(processorCount);

  for i = 1:processorCount
    for j = (i + 1):processorCount
      C(i, j) = exp(-abs(distance(i) - distance(j)) / correlationLength);
      C(j, i) = C(i, j);
    end
  end

  mapping = decomposeSVD(C, threshold);

  %
  % Take into account the grobal parameter.
  %
  mapping = [ sqrt(1 - globalPortion) * mapping, ...
    sqrt(globalPortion) * ones(processorCount, 1) ];
end

function mapping = decomposeSVD(C, threshold)
  [ V, L ] = pcacov(C);

  [ ~, L, I ] = Utils.chooseSignificant(L, threshold);
  V = V(:, I);

  mapping = V * diag(sqrt(L));
end
