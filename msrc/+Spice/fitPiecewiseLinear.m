function [ T, B ] = fitPiecewiseLinear(f, Lnom, Tshift)
  if nargin < 3, Tshift = 0; end

  TT = Utils.toKelvin([ 0 20 40 60 80 100 120 140 160 180 200 ])';

  count = length(TT) - 1;

  T = zeros(count, 2);
  B = zeros(count, 2);

  for i = 1:count
    T0 = (TT(i):0.1:TT(i + 1))';

    Y = f(Lnom, T0);
    X = [ ones(size(T0)) (T0 - Tshift) ];

    T(i, :) = [ TT(i) TT(i + 1) ] - Tshift;
    B(i, :) = (mldivide(X' * X, X' * Y))';
  end
end
