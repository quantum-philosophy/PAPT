init;

methodSet = { 'GaussHermiteProbabilists', 'GaussHermitePhysicists', 'KronrodPatterson' };
orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
coreSet = [ 2 4 8 16 32 ];

methodCount = length(methodSet);
coreCount = length(coreSet);
orderCount = length(orderSet);

sdimSet = zeros(coreCount, 1);

for i = 1:coreCount
  floorplan = Utils.resolveTest(coreSet(i));
  C = PrincipalComponent.perform(floorplan);
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

  fprintf('%15s', 'Order');
  for i = 1:coreCount
    fprintf('%15s', sprintf('Cores %d(%d)', coreSet(i), sdimSet(i)));
  end
  fprintf('\n');

  for i = 1:orderCount
    polynomialOrder = orderSet(i);

    fprintf('%15d', polynomialOrder);

    for j = 1:coreCount
      sdim = sdimSet(j);

      quadratureOrder = ...
        Quadrature.(method).polynomialOrderToQuadratureOrder(polynomialOrder);

      if sparse
        points = ...
          Quadrature.(method).countSparseGridPoints(sdim, quadratureOrder);
      else
        points = ...
          Quadrature.(method).countTensorProductPoints(sdim, quadratureOrder);
      end

      fprintf('%15d', points);
    end
    fprintf('\n');
  end
  fprintf('\n');
end
