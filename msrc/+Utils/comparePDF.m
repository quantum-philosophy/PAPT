function error = comparePDF(outMC, outPC, time, labels)
  [ samples, ddim, tdim ] = size(outMC);

  points = 100;
  draw = false;

  if nargin > 2
    draw = true;

    if nargin < 4
      labels = {};
      for i = 1:ddim
        labels{end + 1} = sprintf('PE %d', i);
      end
    end

    figure;
    title('Probability Density Error');
  end

  error = zeros(ddim, tdim);

  done = 0;
  total = ddim * tdim;

  h = ibar('Comparison of PDFs: step %d out of %d.', total);

  for i = 1:ddim
    for j = 1:tdim
      mc = outMC(:, i, j);
      pc = outPC(:, i, j);

      minMC = min(mc);
      maxMC = max(mc);

      minPC = min(pc);
      maxPC = max(pc);

      x = linspace(min(minMC, minPC), max(maxMC, maxPC), points);

      densityMC = histc(mc, x) / samples;
      densityPC = histc(pc, x) / samples;

      % densityMC = ksdensity(mc, x);
      % densityPC = ksdensity(pc, x);

      error(i, j) = Utils.NRMSE(densityMC, densityPC);

      done = done + 1;
      increase(h);
    end

    if draw
      color = Utils.pickColor(i);
      line(time, error(i, :) * 100, 'Color', color);
    end
  end

  if draw
    legend(labels{:});
    xlabel('Time, s');
    ylabel('NRMSE, %');
  end
end
