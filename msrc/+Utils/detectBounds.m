function [ left, right ] = detectBounds(mc, pc);
  minMC = min(mc);
  minPC = min(pc);

  maxMC = max(mc);
  maxPC = max(pc);

  expMC = mean(mc);
  expPC = mean(pc);

  stdMC = sqrt(var(mc));
  stdPC = sqrt(var(pc));

  leftMC = max(minMC, expMC - 3 * stdMC);
  leftPC = max(minPC, expPC - 3 * stdPC);

  rightMC = min(maxMC, expMC + 3 * stdMC);
  rightPC = min(maxPC, expPC + 3 * stdPC);

  left = min(leftMC, leftPC);
  right = max(rightMC, rightPC);
end
