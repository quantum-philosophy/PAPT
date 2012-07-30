function x = constructLinearSpace(mc, pc)
  [ left, right ] = Utils.detectBounds(mc, pc);
  points = max((right - left) / 0.1, 100);
  x = linspace(left, right, points);
end
