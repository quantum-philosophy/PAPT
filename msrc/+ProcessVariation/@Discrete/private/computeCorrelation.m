function C = computeCorrelation(floorplan)
  D = dlmread(floorplan, '', 0, 1);

  W = D(:, 1);
  H = D(:, 2);
  X = D(:, 3);
  Y = D(:, 4);

  cores = size(D, 1);

  dieW = max(X + W);
  dieH = max(Y + H);

  dieX = dieW / 2;
  dieY = dieH / 2;

  coreX = X + W / 2;
  coreY = Y + H / 2;

  corrLength = max(dieW, dieH) / 2;

  distance = sqrt((dieX - coreX).^2 + (dieY - coreY).^2);

  C = eye(cores);

  for i = 1:cores
    for j = (i + 1):cores
      C(i, j) = exp(-abs(distance(i) - distance(j)) / corrLength);
      C(j, i) = C(i, j);
    end
  end
end
