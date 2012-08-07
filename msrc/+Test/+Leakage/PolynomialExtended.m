init;

f = Spice.fitExponentPolynomial('inverter_45nm', ...
  [ 1 2 ], [ 1, 0.7, 0; 1, 1, 1 ], true);

Lnom = HotSpot.Base.Lnom;
Ldev = PrincipalComponent.deviationRatio * Lnom;

Tref = Utils.toKelvin(27);

T = linspace(Tref, Tref + 123, 100);
L = linspace(Lnom - 4 * Ldev, Lnom + 4 * Ldev, 100);

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

i = f(l, t) / f(HotSpot.Base.Lnom, Tref);

hold on;

plot3(l, t, i, ...
  'LineStyle', 'None', ...
  'Marker', 'o', ...
  'MarkerEdgeColor', 'none', ...
  'MarkerSize', 2, ...
  'MarkerFaceColor', 'g');

grid on;
