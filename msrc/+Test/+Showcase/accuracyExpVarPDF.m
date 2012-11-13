setup;

fontName = 'CMU_Sans_Serif';

data = [
  1  8.64  0.38  35.97;
  2  6.14  0.13  13.96;
  3  3.76  0.14   5.36;
  4  2.05  0.14   2.25;
  5  1.49  0.14   1.64;
  6  1.39  0.14   2.33;
  % 7  1.55  0.15   3.66;
];

order = data(:, 1);
pdf = data(:, 2);
exp = data(:, 3);
var = data(:, 4);

figure;

line(order, exp, ...
  'Color', Color.pick(1), ...
  'LineWidth', 2, 'Marker', 'v', 'MarkerSize', 8);

line(order, var, ...
  'Color', Color.pick(2), ...
  'LineWidth', 2, 'Marker', 's', 'MarkerSize', 8);

line(order, pdf, ...
  'Color', Color.pick(3), ...
  'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);

set(gca, 'XTick', order);

set(gca, ...
  'FontName', fontName, ...
  'FontSize', 14);

legend('Expectation', 'Variance', 'PDF');

title('Error Measurements', ...
  'FontName', fontName, ...
  'FontSize', 20, ...
  'FontWeight', 'normal');

xlabel('Polynomial order', ...
  'FontName', fontName, ...
  'FontSize', 16);

ylabel('Error, %', ...
  'FontName', fontName, ...
  'FontSize', 16);
