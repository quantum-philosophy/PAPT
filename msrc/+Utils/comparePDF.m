function error = comparePDF(timeLine, outMC, outPC, labels, points)
  [ samples, ddim, tdim ] = size(outMC);

  if nargin < 4
    labels = {};
    for i = 1:ddim
      labels{end + 1} = sprintf('PE %d', i);
    end
  end

  if nargin < 5
    points = 100;
  end

  error = zeros(ddim, tdim);

  figure;
  title('Probability Density Error');

  done = 0;
  total = ddim * tdim;

  h = waitbar(0, sprintf('Comparison of PDFs: %d out of %d.', done, total));

  for i = 1:ddim
    for j = 1:tdim
      mc = outMC(:, i, j);
      pc = outPC(:, i, j);

      minMC = min(mc);
      maxMC = max(mc);

      minPC = min(pc);
      maxPC = max(pc);

      x = linspace(max(minMC, minPC), min(maxMC, maxPC), points);

      densityMC = ksdensity(mc, x);
      densityPC = ksdensity(pc, x);
      error(i, j) = Utils.NRMSE(densityMC, densityPC);

      done = done + 1;
      h = waitbar(done / total, h, ...
        sprintf('Comparison of PDFs: %d out of %d.', done, total));
    end

    color = Utils.pickColor(i);
    line(timeLine, error(i, :), 'Color', color);
  end

  close(h);

  legend(labels{:});
end
