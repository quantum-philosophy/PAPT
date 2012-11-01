setup;

quadratureRule = { 'GaussHermite' };
polynomialOrder = [ 1 2 3 4 5 ];
processorCount  = [ 2 4 8 16 32 ];

dimension = zeros(length(processorCount), 1);

for i = 1:length(processorCount)
  floorplan = File.join( ...
    File.trace, '..', 'Assets', sprintf('%02d.flp', processorCount(i)));
  [ ~, dimension(i) ] = ProcessVariation.analyze( ...
    floorplan, LeakagePower.Lnom);
end

for k = 1:(2 * length(quadratureRule))
  ruleName = quadratureRule{ceil(k / 2)};
  sparse = mod(k, 2) == 1;

  if sparse
    fprintf('Quadrature rule: %s (SG)\n', ruleName);
    method = 'sparse';
  else
    fprintf('Quadrature rule: %s (TP)\n', ruleName);
    method = 'tensor';
  end

  fprintf('%10s%10s', 'PC order', 'QD order');
  for i = 1:length(processorCount)
    fprintf('%15s', sprintf('%d / %d', processorCount(i), dimension(i)));
  end
  fprintf('\n');

  for i = 1:length(polynomialOrder)
    order = polynomialOrder(i) + 1;

    fprintf('%10d%10d', polynomialOrder(i), order);

    for j = 1:length(processorCount)
      quadrature = Quadrature('dimension', dimension(j), ...
        'order', order, 'ruleName', ruleName, 'method', method);
      fprintf('%15d', quadrature.nodeCount);
    end
    fprintf('\n');
  end
  fprintf('\n');
end
