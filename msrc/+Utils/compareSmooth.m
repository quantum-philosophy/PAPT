function compareSmooth(outMC, outPC, labels)
  [ samples, ddim ] = size(outMC);

  figure;

  for i = 1:ddim
    subplot(1, ddim, i);
    title('Probability Density');

    mc = outMC(:, i);
    pc = outPC(:, i);

    x = Utils.constructLinearSpace(mc, pc);

    density1 = ksdensity(outMC(:, i), x);
    line(x, density1, 'Color', Utils.pickColor(1));

    density2 = ksdensity(outPC(:, i), x);
    line(x, density2, 'Color', Utils.pickColor(2));

    xlim([ min(x), max(x) ]);

    e = Utils.NRMSE(density1, density2) * 100;
    legend(labels{:});
  end
end
