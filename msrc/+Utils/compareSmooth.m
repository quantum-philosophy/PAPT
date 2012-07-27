function compareSmooth(x, out1, out2, labels)
  [ samples, ddim ] = size(out1);

  if size(out2, 1) ~= samples || size(out2, 2) ~= ddim
    error('The dimensions do not match each other.');
  end

  figure;

  for i = 1:ddim
    subplot(1, ddim, i);

    density1 = ksdensity(out1(:, i), x);
    line(x, density1, 'Color', Utils.pickColor(1));

    density2 = ksdensity(out2(:, i), x);
    line(x, density2, 'Color', Utils.pickColor(2));

    xlim([ min(x), max(x) ]);

    e = Utils.NRMSE(density1, density2);
    title(sprintf('Probability Density (NRMSE %.2e)', e));
    legend(labels{:});
  end
end
