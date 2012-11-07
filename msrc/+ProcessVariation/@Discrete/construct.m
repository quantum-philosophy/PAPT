function C = construct(this, floorplan, options)
  D = dlmread(floorplan, '', 0, 1);

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

  C = performPCA(C, ...
    options.get('reduction', 'adjustable'), ...
    options.get('threshold', this.threshold));
end

function P = performPCA(M, reduction, threshold)
  [ P, L, E ] = pcacov(M);

  switch lower(reduction)
  case 'adjustable'
    keep = min(find((cumsum(E) / sum(E) - threshold) > 0));
    if isempty(keep), keep = size(P, 1); end
  case 'none'
    keep = size(P, 1);
  otherwise
    error('The reduction method is unknown.');
  end

  P = P(:, 1:keep) * diag(sqrt(L(1:keep)));
end
