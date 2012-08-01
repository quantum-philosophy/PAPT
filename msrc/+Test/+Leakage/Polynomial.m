init;

Spice.fitExponentPolynomial('inverter_45nm', ...
  [ 1 2 ], [ 1, 0.7, 0; 1, 1, 1 ], true);

Spice.fitPolynomial('inverter_45nm', ...
  [ 3 2 ], true);
