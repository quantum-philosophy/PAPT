function f = fitExponentPolynomial(name, order, scale, draw)
  if nargin < 2 || numel(order) == 0, order = [ 2 2 ]; end

  assert(numel(order) == 2, 'The order vector is invalid.');

  maxOrder = max(order);

  if nargin < 3 || numel(scale) == 0, scale = ones(2, maxOrder + 1); end

  assert(size(scale, 1) == 2 && size(scale, 2) == maxOrder + 1, ...
    'The scaling matrix is invalid.');

  if nargin < 4, draw = false; end

  %
  % Compute a stamp for the order vector.
  %
  o = sprintf('%d %d', order(1), order(2));

  %
  % Compute a stamp for the scaling matrix.
  %
  s = '';
  for k = 1:2
    if k > 1, s = [ s, ' ' ]; end
    for i = 1:(maxOrder + 1)
      if i > 1, s = [ s, ' ' ]; end
      s = [ s, sprintf('%.2f', scale(k, i)) ];
    end
  end

  filename = [ 'Leakage(ExponentPolynomial_', name, ')', ...
    '_o(', regexprep(o, ' ', '_'), ')', ...
    '_s(', regexprep(s, ' ', '_'), ').mat' ];

  filename = Utils.resolvePath(filename, 'cache');

  if exist(filename, 'file')
    load(filename);
  else
    debug({ 'Construction of a new SPICE fit.' }, ...
          { '  Circuit name: %s', name }, ...
          { '  Type: exponent of a polynomial' }, ...
          { '  Order of L and T: %s', o }, ...
          { '  Scale of L and T: %s', s });

    D = dlmread(Utils.resolvePath([ name, '.leak' ], 'test'), '\t', 1, 0);
    Ldata = D(:, 1);
    Tdata = Utils.toKelvin(D(:, 2));
    Idata = D(:, 3);

    [ fitresult, gof, LTmean, LTstd ] = ...
      Spice.doPolynomialFit(Ldata, Tdata, log(Idata), order);

    cvals = coeffvalues(fitresult);
    cnames = coeffnames(fitresult);

    Lsym = sympoly('L');
    Tsym = sympoly('T');

    Lnorm = (Lsym - LTmean(1)) / LTstd(1);
    Tnorm = (Tsym - LTmean(2)) / LTstd(2);

    logI = sympoly(0);

    for i = 1:numel(cnames)
      attrs = regexp(cnames{i}, '^p(\d)(\d)$', 'tokens');

      Lorder = str2num(attrs{1}{1});
      Torder = str2num(attrs{1}{2});

      s = scale(1, Lorder + 1) * scale(2, Torder + 1);

      logI = logI + s * cvals(i) * Lnorm^Lorder * Tnorm^Torder;
    end

    [ arguments, body ] = Utils.toFunctionString(logI, Lsym, Tsym);
    f =  str2func([ '@(', arguments, ')exp(', body, ')' ]);

    save(filename, 'f', 'Ldata', 'Tdata', 'Idata');
  end

  if ~draw, return; end

  norm = f(HotSpot.Base.Lnom, Utils.toKelvin(27));

  Idata = Idata / norm;

  figure;
  h = subplot(1, 1, 1);

  Luni = sort(unique(Ldata));
  Tuni = sort(unique(Tdata));

  [ L, T ] = meshgrid(Luni, Tuni);

  I = griddata(Ldata, Tdata, Idata, L, T);

  mesh(L, T, I);

  line(Ldata, Tdata, f(Ldata, Tdata) / norm, ...
      'LineStyle', 'None', ...
      'Marker', 'o', ...
      'MarkerEdgeColor', 'w', ...
      'MarkerFaceColor', 'b', ...
      'Parent', h);

  title('Exponent of Polynomial Fit');
  xlabel('L');
  ylabel('T');
  zlabel('I');

  grid on;
  view(10, 10);
end
