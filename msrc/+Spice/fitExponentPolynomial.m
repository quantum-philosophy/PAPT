function f = fitExponentPolynomial(name, order, scale, draw)
  if nargin < 2 || numel(order) == 0, order = [ 2 2 ]; end
  if nargin < 3 || numel(scale) == 0, scale = 1.0; end
  if nargin < 4, draw = false; end

  filename = [ 'Leakage_exponent_polynomial_', name, ...
    '_l', num2str(order(1)), '_t', num2str(order(2)), ...
    sprintf('_s%.2f', scale), '.mat' ];
  filename = Utils.resolvePath(filename, 'cache');

  if exist(filename, 'file')
    load(filename);
  else
    debug({ 'Construction of a new SPICE fit.' }, ...
          { '  Circuit name: %s', name }, ...
          { '  Type: exponent of a polynomial' }, ...
          { '  Order of L: %d', order(1) }, ...
          { '  Order of T: %d', order(2) }, ...
          { '  Scale: %.2f', scale });

    D = dlmread(Utils.resolvePath([ name, '.leak' ], 'test'), '\t', 1, 0);
    Ldata = D(:, 1);
    Tdata = Utils.toKelvin(D(:, 2));
    Idata = D(:, 3);

    [ fitresult, gof, LTmean, LTstd ] = ...
      Spice.doPolynomialFit(Ldata, Tdata, log(Idata), order);

    cvals = coeffvalues(fitresult);
    cnames = coeffnames(fitresult);

    Lsym = ipoly('L');
    Tsym = ipoly('T');

    Lnorm = (Lsym - LTmean(1)) / LTstd(1);
    Tnorm = (Tsym - LTmean(2)) / LTstd(2);

    logI = ipoly(0);

    for i = 1:numel(cnames)
      attrs = regexp(cnames{i}, '^p(\d)(\d)$', 'tokens');

      Lorder = str2num(attrs{1}{1});
      Torder = str2num(attrs{1}{2});

      if Lorder > 0 || Torder > 0
        logI = logI + scale * cvals(i) * Lnorm^Lorder * Tnorm^Torder;
      else
        logI = logI + cvals(i) * Lnorm^Lorder * Tnorm^Torder;
      end
    end

    [ arguments, body ] = Utils.toFunctionString(logI, Lsym, Tsym);
    f =  str2func([ '@(', arguments, ')exp(', body, ')' ]);

    save(filename, 'f', 'Ldata', 'Tdata', 'Idata');
  end

  if ~draw, return; end

  h = subplot(1, 1, 1);

  Luni = sort(unique(Ldata));
  Tuni = sort(unique(Tdata));

  [ L, T ] = meshgrid(Luni, Tuni);

  I = griddata(Ldata, Tdata, Idata, L, T);

  mesh(L, T, I);

  line(Ldata, Tdata, f(Ldata, Tdata), ...
      'LineStyle', 'None', ...
      'Marker', 'o', ...
      'MarkerEdgeColor', 'w', ...
      'MarkerFaceColor', 'b', ...
      'Parent', h);

  title('Fitted Surface');
  xlabel('L');
  ylabel('T');
  zlabel('I');

  grid on;
  view(10, 10);
end
