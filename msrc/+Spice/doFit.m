function [ fitresult, gof, L, T, I, mean, std ] = doFit(filename, draw)
  if nargin < 1, filename = Utils.resolvePath('inverter_LT.out'); end
  if nargin < 2, draw = false; end

  M = dlmread(filename, '\t', 1, 0);

  L = M(:, 1);
  T = Utils.toKelvin(M(:, 2));
  I = M(:, 3);

  Tonly = length(unique(L)) == 1;
  Lonly = length(unique(T)) == 1;

  if Lonly
    ft = fittype('poly3');
    opts = fitoptions(ft);
    [ x, y ] = prepareCurveData(L, I);
    vars = 4;
  elseif Tonly
    ft = fittype('poly3');
    opts = fitoptions(ft);
    [ x, y ] = prepareCurveData(T, I);
    vars = 4;
  else
    ft = fittype('poly33');
    opts = fitoptions(ft);
    [ x, y, z ] = prepareSurfaceData(L, T, I);
    x = [ x, y ];
    y = z;
    vars = 10;
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
