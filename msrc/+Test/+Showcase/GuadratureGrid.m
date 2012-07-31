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
requiredOrder = order + 1;

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
fprintf('Integration order:   %d\n', requiredOrder);

fprintf('\n');

[ nodes, weights, points ] = ...
  GaussianQuadrature.Physicists.constructSparseGrid(sdim, requiredOrder);
fprintf('Sparse physicists:\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  GaussianQuadrature.Probabilists.constructSparseGrid(sdim, requiredOrder);
fprintf('Sparse probabilists:\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  GaussianQuadrature.Physicists.constructTensorProduct(sdim, requiredOrder);
fprintf('Tensor physicists:\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');

[ nodes, weights, points ] = ...
  GaussianQuadrature.Probabilists.constructTensorProduct(sdim, requiredOrder);
fprintf('Tensor probabilists:\n');
fprintf('  Total points:      %d\n', points);
fprintf('  Negative points:   %d\n', length(find(weights < 0)));
fprintf('  Integration:       %e\n', sum(f(nodes).^2 .* weights));
fprintf('\n');
