function plotExpStd(time, exp, var)
  figure;

  std = sqrt(var);
  count = size(exp, 1);

  for i = 1:count
    color = Utils.pickColor(i);
    line(time, exp(i, :), 'Color', color);
    line(time, exp(i, :) + std(i, :), 'Color', color, 'LineStyle', '--');
    line(time, exp(i, :) - std(i, :), 'Color', color, 'LineStyle', '--');
  end

  xlabel('Time, s');
  ylabel('Temperature, C');
end
