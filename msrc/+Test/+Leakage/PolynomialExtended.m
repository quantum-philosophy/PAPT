init;

f = Spice.fitExponentPolynomial('inverter_45nm', ...
  [ 1 2 ], [ 1, 0.7, 0; 1, 1, 1 ], true);

Lnom = 45e-9;
Ldev = 0.05 * Lnom;

T = Utils.toKelvin(linspace(27, 150, 100));
L = linspace(Lnom - 3 * Ldev, Lnom + 3 * Ldev, 100);

l = zeros(100 * 100, 1);
t = zeros(100 * 100, 1);

k = 0;
for i = 1:100
  for j = 1:100
    k = k + 1;
    l(k) = L(i);
    t(k) = T(j);
  end
end

i = f(l, t);

hold on;

plot3(l, t, i, ...
  'LineStyle', 'None', ...
  'Marker', 'o', ...
  'MarkerEdgeColor', 'none', ...
  'MarkerSize', 2, ...
  'MarkerFaceColor', 'g');

grid on;
