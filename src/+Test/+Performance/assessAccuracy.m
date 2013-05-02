function assessAccuracy
  setup;
  rng(1);

  stepCount = 1e2;
  chaosSampleCount = 1e5;

  comparisonOptions = Options('method', 'histogram', ...
    'range', 'unbounded', 'function', 'pdf');

  display(comparisonOptions, 'Comparison options');

  orderSet       = [ 1 2 3 4 5 6 7 ];
  sampleCountSet = [ 1e2 1e3 1e4 1e5 ];

  pick = [ 7 1e4 50 ];

  orderCount = length(orderSet);
  sampleCount = length(sampleCountSet);

  options = Test.configure('stepCount', stepCount);

  %
  % Monte Carlo simulation.
  %
  mc = Temperature.MonteCarlo.Transient(options.temperatureOptions);

  [ ~, output ] = mc.compute(options.dynamicPower, ...
    'sampleCount', max(sampleCountSet), 'verbose', true);

  mcTDATA = output.Tdata(randperm(max(sampleCountSet)), :, :);

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

    chaos = Temperature.Chaos.Transient(options.temperatureOptions, ...
      options.chaosOptions);

    [ Texp, output ] = chaos.compute(options.dynamicPower);
    Tdata = chaos.sample(output.coefficients, chaosSampleCount);

    for j = 1:sampleCount
      errorExp(i, j) = Error.computeRMSE(mcTexp{j}, Texp);
      errorVar(i, j) = Error.computeRMSE(mcTvar{j}, output.Tvar);
      errorPDF(i, j) = Data.compare(mcTdata{j}, Tdata, ...
        comparisonOptions);

      if orderSet(i) == pick(1) && sampleCountSet(j) == pick(2)
        Data.compare(mcTdata{j}(:, :, pick(3)), Tdata(:, :, pick(3)), ...
          comparisonOptions, 'draw', true, 'layout', 'separate');
      end

      fprintf('%15.6f', errorPDF(i, j));
    end

    fprintf(' | ');

    for j = 1:sampleCount
      fprintf('%15.6f', errorExp(i, j));
    end

    fprintf(' | ');

    for j = 1:sampleCount
      fprintf('%15.6f', errorVar(i, j));
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
