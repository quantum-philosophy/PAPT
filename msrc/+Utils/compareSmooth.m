function compareSmooth(outMC, outPC, options)
  [ samples, ddim ] = size(outMC);

  figure;

  for i = 1:ddim
    subplot(1, ddim, i);
    title('Empirical Probability Density');

    mc = outMC(:, i);
    pc = outPC(:, i);

    if isfield(options, 'x')
      x = options.x;
    else
      x = Utils.constructLinearSpace(mc, pc, options);
    end

    density1 = ksdensity(mc, x);
    line(x, density1, 'Color', 'b');

    density2 = ksdensity(pc, x);
    line(x, density2, 'Color', 'g');

    if isfield(options, 'exact')
      f = options.exact;
      density0 = f(x, i);
      line(x, density0, 'Color', 'r');
    end

    legend(options.labels{:});
  end
end
