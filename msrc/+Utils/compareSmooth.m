function compareSmooth(outMC, outPC, labels)
  [ samples, ddim ] = size(outMC);

  if size(outPC, 1) ~= samples || size(outPC, 2) ~= ddim
    error('The dimensions do not match each other.');
  end

  points = 100;

  figure;

  for i = 1:ddim
    subplot(1, ddim, i);

    mc = outMC(:, i);
    pc = outPC(:, i);

    minMC = min(mc);
    maxMC = max(mc);

    minPC = min(pc);
    maxPC = max(pc);

    x = linspace(max(minMC, minPC), min(maxMC, maxPC), points);

    density1 = ksdensity(outMC(:, i), x);
    line(x, density1, 'Color', Utils.pickColor(1));

    density2 = ksdensity(outPC(:, i), x);
    line(x, density2, 'Color', Utils.pickColor(2));

    xlim([ min(x), max(x) ]);

    e = Utils.NRMSE(density1, density2) * 100;
    title(sprintf('Probability Density (NRMSE %.2f%%)', e));
    legend(labels{:});
  end
end
