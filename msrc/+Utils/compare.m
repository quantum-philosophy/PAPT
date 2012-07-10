function compare(out1, out2, labels, f)
  exact = (nargin > 3);

  [ samples, ddim ] = size(out1);

  if size(out2, 1) ~= samples || size(out2, 2) ~= ddim
    error('The dimensions do not match each other.');
  end

  min1 = min(out1);
  max1 = max(out1);

  min2 = min(out2);
  max2 = max(out2);

  figure;

  rows = 2;
  cols = ddim;

  for i = 1:ddim
    subplot(rows, cols, i);

    x = linspace(max(min1(i), min2(i)), min(max1(i), max2(i)), 200);

    if exact
      density0 = f(x, i);
      line(x, density0, 'Color', 'r');
    end

    density1 = ksdensity(out1(:, i), x);
    line(x, density1, 'Color', 'b');

    density2 = ksdensity(out2(:, i), x);
    line(x, density2, 'Color', 'g');

    xlim([ x(1), x(end) ]);

    title('PDF');
    if exact
      legend('Exact', labels{:});
    else
      legend(labels{:});
    end

    subplot(rows, cols, cols + i);

    n = length(x);

    if exact
      line(x, density0 - density1, 'Color', 'g');
      line(x, density0 - density2, 'Color', 'b');

      rmse1 = sqrt(sum((density0 - density1).^2) / n);
      rmse2 = sqrt(sum((density0 - density2).^2) / n);

      title('RMSE');
      legend(sprintf('%s (%.4e)', labels{1}, rmse1), ...
        sprintf('%s (%.4e)', labels{2}, rmse2));
    else
      line(x, density1 - density2, 'Color', 'r');

      rmse = sqrt(sum((density1 - density2).^2) / n);

      title(sprintf('RMSE %.4e', rmse));
      legend([ labels{1}, ' - ', labels{2} ]);
    end

    xlim([ x(1), x(end) ]);
  end
end
