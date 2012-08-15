init;

methodSet = { ...
  'GaussHermiteProbabilists', ...
  'GaussHermitePhysicists', ...
  'GaussJacobi', ...
  'KronrodPatterson' };

orderSet = [ 1 2 3 4 5 ];
coreSet = [ 2 4 8 16 32 ];

methodCount = length(methodSet);
coreCount = length(coreSet);
orderCount = length(orderSet);

sdimSet = zeros(coreCount, 1);

for i = 1:coreCount
  floorplan = Utils.resolveTest(coreSet(i));
  C = PrincipalComponent.perform(floorplan, HotSpot.Base.Lnom);
  sdimSet(i) = size(C, 2);
end

for k = 1:(2 * methodCount)
  method = methodSet{ceil(k / 2)};

  sparse = mod(k, 2) == 1;

  if sparse
    fprintf('Method: %s (SG)\n', method);
  else
    fprintf('Method: %s (TP)\n', method);
  end

  fprintf('%10s%10s%10s', 'PC Order', 'QD Order', 'QD Level');
  for i = 1:coreCount
    fprintf('%15s', sprintf('Cores %d(%d)', coreSet(i), sdimSet(i)));
  end
  fprintf('\n');

  for i = 1:orderCount
    chaosOrder = orderSet(i);
    quadratureOrder = chaosOrder + 1;
    quadratureLevel = chaosOrder; % ceil(log2(quadratureOrder + 1) - 1);

    fprintf('%10d%10d%10d', chaosOrder, quadratureOrder, quadratureLevel);

    for j = 1:coreCount
      sdim = sdimSet(j);

      if sparse
        points = Quadrature.(method).countSparseGridPoints(...
          sdim, quadratureOrder, quadratureLevel);
      else
        points = Quadrature.(method).countTensorProductPoints(...
          sdim, quadratureOrder);
      end

      fprintf('%15d', points);
    end
    fprintf('\n');
  end
  fprintf('\n');
end
