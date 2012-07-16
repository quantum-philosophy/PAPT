function [ fitresult, gof ] = fit(filename)
  if nargin < 1, filename = Utils.resolvePath('nmos.out'); end

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
    vars = 10;
  end

  opts.Normalize = 'on';
  opts.Lower = ones(1, vars) * -Inf;
  opts.Upper = ones(1, vars) * Inf;

  if Lonly
    [ fitresult, gof ] = fit(x, y, ft, opts);
    plot(fitresult, x, y);
    xlabel('L');
    ylabel('I');
    legend('Original', 'Fitted');
  elseif Tonly
    [ fitresult, gof ] = fit(x, y, ft, opts);
    plot(fitresult, x, y);
    xlabel('T');
    ylabel('I');
    legend('Original', 'Fitted');
  else
    [ fitresult, gof ] = fit([ x, y ], z, ft, opts);
    plot(fitresult, [ x, y ], z);
    xlabel('L');
    ylabel('T');
    zlabel('I');
    legend('Fitted', 'Original');
  end
  grid on;
end
