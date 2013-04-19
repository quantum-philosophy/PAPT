function assessAccuracy
  setup;
  rng(1);

  stepCount = 1e2;
  chaosSampleCount = 1e5;

  comparisonOptions = Options( ...
    'method', 'histogram', 'range', 'unbounded', 'function', 'pdf');

  display(comparisonOptions, 'Comparison options');

  orderSet       = [ 1 2 3 4 5 6 7 ];
  sampleCountSet = [ 1e2 1e3 1e4 1e5 ];

  pick = [ 0 0 ];

  orderCount = length(orderSet);
  sampleCount = length(sampleCountSet);

  options = Test.configure('stepCount', stepCount);

  %
  % Monte Carlo simulation.
  %
  mc = HotSpot.MonteCarlo(options.hotspotOptions, ...
    'sampleCount', max(sampleCountSet), 'verbose', true);

  [ ~, ~, mcTDATA ] = mc.computeInParallel( ...
    options.powerProfile, options.leakage);

  mcTDATA = mcTDATA(randperm(max(sampleCountSet)), :, :);

  mcTexp = cell(sampleCount, 1);
  mcTvar = cell(sampleCount, 1);
  mcTdata = cell(sampleCount, 1);

  for i = 1:sampleCount
    mcTdata{i} = mcTDATA(1:sampleCountSet(i), :, :);
    mcTexp{i} = mean(mcTdata{i}, 1);
    mcTvar{i} = var(mcTdata{i}, [], 1);
  end

  errorExp = zeros(orderCount, sampleCount);
  errorVar = zeros(orderCount, sampleCount);
  errorPDF = zeros(orderCount, sampleCount);

  printHeader(sampleCountSet);

  %
  % Polynomial chaos expansion.
  %
  for i = 1:orderCount
    options.chaosOptions.order = orderSet(i);
    options.chaosOptions.quadratureOptions.polynomialOrder = orderSet(i);

    fprintf('%5d | ', orderSet(i));

    chaos = HotSpot.Chaos(options.hotspotOptions, options.chaosOptions);

    [ Texp, Tvar, coefficients ] = chaos.compute( ...
      options.powerProfile, options.leakage);

    Tdata = chaos.sample(coefficients, chaosSampleCount);

    for j = 1:sampleCount
      errorExp(i, j) = Error.computeNRMSE(mcTexp{j}, Texp) * 100;
      errorVar(i, j) = Error.computeNRMSE(mcTvar{j}, Tvar) * 100;

      if all([ i, j ] == pick)
        errorPDF(i, j) = Data.compare(mcTdata{j}, Tdata, ...
          comparisonOptions, 'draw', true) * 100;
      else
        errorPDF(i, j) = Data.compare(mcTdata{j}, Tdata, ...
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
