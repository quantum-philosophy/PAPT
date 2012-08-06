function x = ibetarnd(alpha, beta, a, b, rows, cols)
  x = (b - a) * betarnd(alpha, beta, rows, cols) + a;
end
