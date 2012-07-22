function M = replicate(T, steps)
  [ rows, cols ] = size(T);
  M = zeros(rows, steps);
  packed = 0;
  while (packed < steps)
    topack = min(steps - packed, cols);
    M(:, (packed + 1):(packed + topack)) = T(:, 1:topack);
    packed = packed + topack;
  end
end
