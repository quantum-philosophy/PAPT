function [ fitresult, gof, L, T, I, mean, std ] = doFit(filename, order, draw)
  if nargin < 2, order = [ 3 3 ]; end
  if nargin < 3, draw = false; end

  M = dlmread(filename, '\t', 1, 0);

  L = M(:, 1);
  T = Utils.toKelvin(M(:, 2));
  I = M(:, 3);

  Tonly = length(unique(L)) == 1;
  Lonly = length(unique(T)) == 1;

  warn(filename, order, Lonly, Tonly);

  if Lonly
    order = order(1);
    ft = fittype(sprintf('poly%d', order));
    opts = fitoptions(ft);
    [ x, y ] = prepareCurveData(L, I);
    vars = order + 1;
  elseif Tonly
    order = order(1);
    ft = fittype(sprintf('poly%d', order));
    opts = fitoptions(ft);
    [ x, y ] = prepareCurveData(T, I);
    vars = order + 1;
  else
    ft = fittype(sprintf('poly%d%d', order(1), order(2)));
    opts = fitoptions(ft);
    [ x, y, z ] = prepareSurfaceData(L, T, I);
    x = [ x, y ];
    y = z;
    vars = order(1) * order(2) + 1;
  end

  opts.Normalize = 'off';

  [ x, mean, std ] = curvefit.normalize(x);

  opts.Lower = ones(1, vars) * -Inf;
  opts.Upper = ones(1, vars) * Inf;

  [ fitresult, gof ] = fit(x, y, ft, opts);

  if ~draw, return; end

  plot(fitresult, x, y);

  if Lonly
    xlabel('L');
    ylabel('I');
    legend('Original', 'Fitted');
  elseif Tonly
    xlabel('T');
    ylabel('I');
    legend('Original', 'Fitted');
  else
    xlabel('L');
    ylabel('T');
    zlabel('I');
    legend('Fitted', 'Original');
  end
  grid on;
end

function warn(filename, order, Lonly, Tonly)
  if Lonly
    L = order(1);
    T = 0;
  elseif Tonly
    L = 0;
    T = order(1);
  else
    L = order(1);
    T = order(2);
  end

  debug('------------------------------\n');
  debug('Fitting to new SPICE data:\n');
  debug('  File name: %s\n', filename);
  debug('  Order for channel length: %d\n', L);
  debug('  Order for temperature: %d\n', T);
  debug('------------------------------\n');
end
