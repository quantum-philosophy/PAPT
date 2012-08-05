clc;

n = 10^2;
m = 10^5;

M = randi(100, n, m);
V = randi(100, n, 1);

R0 = zeros(n, m);
R1 = zeros(n, m);

tic
for i = 1:m
  R0(:, i) = M(:, i) .* V;
end
t = toc;
fprintf('for:    %10.6f s\n', t);

tic
R1 = M .* irep(V, 1, m);
t = toc;
fprintf('irep:   %10.6f s (%.2e)\n', t, mean(mean(R1 - R0)));

tic
R1 = M .* repmat(V, 1, m);
t = toc;
fprintf('repmat: %10.6f s (%.2e)\n', t, mean(mean(R1 - R0)));

tic
R1 = bsxfun(@times, M, V);
t = toc;
fprintf('bsxfun: %10.6f s (%.2e)\n', t, mean(mean(R1 - R0)));