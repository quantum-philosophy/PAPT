function f = fit(name, order, draw)
  if nargin < 2 || numel(order) == 0, order = [ 3 3 ]; end
  if nargin < 3, draw = false; end

  filename = [ 'LT_', name, ...
    '_l', num2str(order(1)), '_t', num2str(order(2)), '.mat' ];
  filename = Utils.resolvePath(filename, 'cache');

  if exist(filename, 'file')
    load(filename);
  else
    [ fitresult, gof, Ldata, Tdata, Idata, LTmean, LTstd ] = ...
      Spice.doFit(Utils.resolvePath([ name, '.leak' ]), order);

    cvals = coeffvalues(fitresult);
    cnames = coeffnames(fitresult);

    Lsym = ipoly('L');
    Tsym = ipoly('T');

    Lnorm = (Lsym - LTmean(1)) / LTstd(1);
    Tnorm = (Tsym - LTmean(2)) / LTstd(2);

    I = ipoly(0);

    for i = 1:numel(cnames)
      attrs = regexp(cnames{i}, '^p(\d)(\d)$', 'tokens');

      Lorder = str2num(attrs{1}{1});
      Torder = str2num(attrs{1}{2});

      I = I + cvals(i) * Lnorm^Lorder * Tnorm^Torder;
    end

    f = Utils.toFunction(I, Lsym, Tsym);

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
end
