function compareHistogram(x, out1, out2, labels)
  [ samples, ddim ] = size(out1);

  if size(out2, 1) ~= samples || size(out2, 2) ~= ddim
    error('The dimensions do not match each other.');
  end

  min1 = min(out1);
  max1 = max(out1);

  min2 = min(out2);
  max2 = max(out2);

  figure;

  for i = 1:ddim
    p = subplot(1, ddim, i);

    hist(p, out1(:, i), x);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75);
    hold on;

    hist(p, out2(:, i), x);
    h1 = findobj(gca, 'Type', 'patch');
    set(h1, 'facealpha', 0.75);

    xlim([ min(x), max(x) ]);

    title('Histogram');
    legend(labels{:});
  end
end
