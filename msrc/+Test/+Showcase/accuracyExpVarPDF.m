setup;

fontName = 'CMU_Sans_Serif';

data = [
  1  9.01  0.65  38.11;
  2  6.24  0.33  17.35;
  3  3.98  0.29   9.12;
  4  2.19  0.27   5.87;
  5  1.53  0.27   4.26;
  6  1.47  0.27   3.22;
  7  1.45  0.26   2.45;
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
