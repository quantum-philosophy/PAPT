setup;

fontName = 'CMU_Sans_Serif';

order = 1:5;

exp = [ 0.16 0.09 0.11 0.11 0.11 ];
var = [ 36.44 11.96 2.19 1.67 1.68 ];
pdf = [ 8.11 6.61 4.11 1.86 1.48 ];

color = [ Color.pick(1); Color.pick(2); Color.pick(3) ];

figure;

line(order, exp, ...
  'Color', color(1, :), ...
  'LineWidth', 2, 'Marker', 'v', 'MarkerSize', 12);

line(order, var, ...
  'Color', color(2, :), ...
  'LineWidth', 2, 'Marker', 's', 'MarkerSize', 12);

line(order, pdf, ...
  'Color', color(3, :), ...
  'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 12);

set(gca, 'XTick', order);

set(gca, ...
  'FontName', fontName, ...
  'FontSize', 14);

legend('Expectation', 'Variance', 'PDF');

xlabel('Polynomial order', ...
  'FontName', fontName, ...
  'FontSize', 16);

ylabel('Error, %', ...
  'FontName', fontName, ...
  'FontSize', 16);
