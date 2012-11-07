function P = performPCA(M, reduction, threshold)
  %
  % Description:
  %
  %   Performs the PCA on the given matrix and shrinks
  %   the number of principal components according to
  %   the threshold, which is constant for now.
  %
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
