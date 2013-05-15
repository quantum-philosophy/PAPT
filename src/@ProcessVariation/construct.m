function construct(this, options)
  kernel = options.kernel;
  threshold = options.threshold;
  globalPortion = options.globalPortion;

  die = options.die;
  D = die.floorplan;

  W = D(:, 1);
  H = D(:, 2);
  X = D(:, 3);
  Y = D(:, 4);

  X = X + W / 2 - die.width / 2;
  Y = Y + H / 2 - die.height / 2;

  I = Utils.constructPairIndex(size(D, 1));
  C = kernel( ...
    [ X(I(:, 1)).'; Y(I(:, 1)).' ], ...
    [ X(I(:, 2)).'; Y(I(:, 2)).' ]);
  C = Utils.symmetrizePairIndex(C, I);

  mapping = decomposeSVD(C, threshold);

  %
  % Take into account the grobal parameter.
  %
  mapping = [ sqrt(1 - globalPortion) * mapping, ...
    sqrt(globalPortion) * ones(die.processorCount, 1) ];

  this.mapping = mapping;
  this.dimensionCount = size(mapping, 2);
end

function mapping = decomposeSVD(C, threshold)
  [ V, L ] = pcacov(C);

  [ ~, L, I ] = Utils.chooseSignificant(L, threshold);
  V = V(:, I);

  mapping = V * diag(sqrt(L));
end
