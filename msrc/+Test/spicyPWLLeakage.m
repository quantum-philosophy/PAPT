init;

Lnom = 45e-9;

f = Spice.fitExponentPolynomial('inverter_45nm', [ 2 2 ], 0.7);
[ T, B ] = Spice.fitPiecewiseLinear(f, Lnom);

TT = [ T(:, 1); T(end, 2) ];

figure;

title('Leakage Current');
X = TT(1):0.1:TT(end);
line(X, f(Lnom, X), 'Color', 'k');
line(TT, f(Lnom, TT), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');

fprintf('%15s%15s\n', 'b0', 'b1');

names = {};

for i = 1:length(T)
  fprintf('%15.4e%15.4e\n', B(i, 1), B(i, 2));
  names{end + 1} = sprintf('[ %.0f, %.0f ]', ...
    Utils.toCelsius(T(i, 1)), Utils.toCelsius(T(i, 2)));

  color = Utils.pickColor(i);
  line(T(i, :), B(i, 1) + B(i, 2) * T(i, :), 'Color', color);
end

legend('Original fit', 'Collocation points', names{:});
