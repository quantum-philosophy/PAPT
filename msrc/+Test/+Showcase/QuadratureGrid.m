init;

%
% The number of stochastic dimensions.
%
sdim = 2;

%
% The maximal total order of polynomials to integrate,
% i.e., the order of the PC expansion.
%
order = 10;

%
% The order (points) of the quadrature rule
% (as if it was one-dimensional).
%
quadratureOrder = order + 1;

%
% The product of two one-dimensional polynomials psi_7 * psi_3.
%
f = @(x) ...
  + 315 .* x(1, :)    .* x(2, :)    ...
  - 105 .* x(1, :)    .* x(2, :).^3 ...
  - 315 .* x(1, :).^3 .* x(2, :)    ...
  + 105 .* x(1, :).^3 .* x(2, :).^3 ...
  +  63 .* x(1, :).^5 .* x(2, :)    ...
  -  21 .* x(1, :).^5 .* x(2, :).^3 ...
  -   3 .* x(1, :).^7 .* x(2, :)    ...
  +        x(1, :).^7 .* x(2, :).^3;

fprintf('Random variables:    %d\n', sdim);
fprintf('Polynomial order:    %d\n', order);
fprintf('Integration order:   %d\n', quadratureOrder);

fprintf('\n');

[ nodes, weights, points ] = ...
  Quadrature.GaussHermitePhysicists.constructSparseGrid(sdim, quadratureOrder);
fprintf('Gauss-Hermite Physicists (SG):\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  Quadrature.GaussHermiteProbabilists.constructSparseGrid(sdim, quadratureOrder);
fprintf('Gauss-Hermite Probabilists (SG):\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  Quadrature.GaussHermitePhysicists.constructTensorProduct(sdim, quadratureOrder);
fprintf('Gauss-Hermite Physicists (TP):\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  Quadrature.GaussHermiteProbabilists.constructTensorProduct(sdim, quadratureOrder);
fprintf('Gauss-Hermite Probabilists (TP):\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  Quadrature.KronrodPatterson.constructSparseGrid(sdim, quadratureOrder);
fprintf('Kronrod-Patterson (SG):\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');
