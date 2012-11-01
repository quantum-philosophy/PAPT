function assessAccuracy
  setup;
  rng(0);

  stepCount = 1e2;
  chaosSampleCount = 1e5;

  comparisonOptions = Options( ...
    'method', 'histogram', 'range', '3sigma', 'function', 'pdf');

  display(comparisonOptions, 'Comparison options');

  orderSet       = [ 1 2 3 4 5 ];
  sampleCountSet = [ 1e2 1e3 1e4 1e5 ];

  pick = [ 0 0 ];

  orderCount = length(orderSet);
  sampleCount = length(sampleCountSet);

  options = configure('stepCount', stepCount);

  printHeader(sampleCountSet);

  %
  % Monte Carlo simulation.
  %
  mcTexp = cell(sampleCount, 1);
  mcTvar = cell(sampleCount, 1);
  mcTdata = cell(sampleCount, 1);

  for i = 1:sampleCount
    numeric = HotSpot.MonteCarlo(options.floorplan, ...
      options.hotspotConfig, options.hotspotLine);

    [ mcTexp{i}, mcTvar{i}, mcTdata{i} ] = ...
      numeric.computeWithLeakageInParallel( ...
      options.powerProfile, options.leakage);

    mcTdata{i} = permute(mcTdata{i}, [ 3, 1, 2 ]);
  end

  errorExp = zeros(orderCount, sampleCount);
  errorVar = zeros(orderCount, sampleCount);
  errorPDF = zeros(orderCount, sampleCount);

  %
  % Polynomial chaos expansion.
  %
  for i = 1:orderCount
    options.chaosOptions.order = orderSet(i);

    fprintf('%5d | ', orderSet(i));

    chaos = HotSpot.StepwiseChaos(options.floorplan, ...
      options.hotspotConfig, options.hotspotLine, options.chaosOptions);

    [ Texp, Tvar, coefficients ] = chaos.computeWithLeakage( ...
      options.powerProfile, options.leakage);

    Tdata = chaos.sample(coefficients, chaosSampleCount);

    for j = 1:sampleCount
      errorExp(i, j) = Error.computeNRMSE(mcTexp{j}, Texp) * 100;
      errorVar(i, j) = Error.computeNRMSE(mcTvar{j}, Tvar) * 100;

      if all([ i, j ] == pick)
        errorPDF(i, j) = compareData(mcTdata{j}, Tdata, ...
          comparisonOptions, 'draw', true) * 100;
      else
        errorPDF(i, j) = compareData(mcTdata{j}, Tdata, ...
          comparisonOptions) * 100;
      end

      fprintf('%15.2f', errorPDF(i, j));
    end

    fprintf(' | ');

    for j = 1:sampleCount
      fprintf('%15.2f', errorExp(i, j));
    end

    fprintf(' | ');

    for j = 1:sampleCount
      fprintf('%15.2f', errorVar(i, j));
    end

    fprintf('\n');
  end
end

function printHeader(sampleCountSet)
  sampleCount = length(sampleCountSet);

  fprintf('\n');

  names = { 'NRMSE(PDF)', 'NRMSE(Exp)', 'NRMSE(Var)' };

  fprintf('%5s | ', '');
  for i = 1:length(names)
    if i > 1, fprintf(' | '); end

    s = names{i};
    start = true;
    while length(s) < 15 * sampleCount
      if start
        s = [ ' ', s ];
        start = false;
      else
        s = [ s, ' ' ];
        start = true;
      end
    end

    fprintf(s);
  end
  fprintf('\n');

  fprintf('%5s | ', 'Order');
  for i = 1:length(names)
    if i > 1, fprintf(' | '); end

    for j = 1:sampleCount
      fprintf('%15s', sprintf('MC %.1e', sampleCountSet(j)));
    end
  end
  fprintf('\n');
end
