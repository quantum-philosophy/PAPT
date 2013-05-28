function plotDensity(filename)
  setup;

  alpha = 0.05;
  mctol = 1e-3;
  shrunkStepCount = 25;
  maximalSampleCount = 1e4;

  load(filename);
  DATA = Tdata;

  [ processorCount, stepCount, sampleCount ] = size(DATA);

  sampleCount = min(sampleCount, maximalSampleCount);
  DATA = Utils.toCelsius(DATA(:, :, 1:sampleCount));

  startStep = 1 + stepCount / 2 - floor(shrunkStepCount / 2);
  stepCount = shrunkStepCount;
  endStep = startStep + stepCount - 1;
  DATA = DATA(:, startStep:endStep, :);

  fprintf('Sample count: %d\n', sampleCount);
  fprintf('Step range:   [ %d, %d ]\n', startStep, endStep);

  k = floor(stepCount / 2);

  figure;

  x = (startStep:endStep)' * 1e-3;

  E = mean(DATA, 3);
  V = var(DATA, [], 3);

  labels = {};

  for i = 1:processorCount
    line(x, E(i, :), ...
      'Color', Color.pick(i), 'LineWidth', 1.5);
    line(x, E(i, :) + sqrt(V(i, :)), ...
      'Color', Color.pick(i), 'LineStyle', '--');
    labels{end + 1} = sprintf('PE %d: expectation', i);
    labels{end + 1} = sprintf('PE %d: deviation', i);
  end

  Plot.label('Time, s', 'Temperature, C');
  Plot.legend(labels{:});

  xlim([ x(1), x(end) ]);

  figure;

  passed1 = 0;
  passed2 = 0;

  labels = {};

  fprintf('%5s %16s %5s %16s\n', 'H1', 'p1', 'H2', 'p2');
  for i = 1:0
    for j = 1:0
      data = squeeze(DATA(j, i, :));

      if i == k
        [ pdf, x ] = ksdensity(data, 'npoints', sampleCount);
        line(x, pdf, 'Color', Color.pick(j), 'LineWidth', 1);
        labels{end + 1} = sprintf('PE %d', j);
      end

      % continue;

      [ tdata, lambda ] = boxcox(data);

      [ H1, p1 ] = jbtest(data, alpha, mctol);
      [ H2, p2 ] = jbtest(tdata, alpha, mctol);

      fprintf('%5d %16.12f %5d %16.12f\n', H1, p1, H2, p2);

      if H1 == 0, passed1 = passed1 + 1; end
      if H2 == 0, passed2 = passed2 + 1; end
    end
  end

  Plot.label('Temperature, C', 'Probability density function');
  Plot.legend(labels{:});

  fprintf('Tests passed for the original data: %d / %d.\n', ...
    passed1, stepCount * processorCount);
  fprintf('Tests passed for the transformed data: %d / %d.\n', ...
    passed2, stepCount * processorCount);

  data = squeeze(DATA(1, k, :));
  tdata = log(data);

  figure;
  subplot(1, 2, 1);
  normplot(data);

  subplot(1, 2, 2);
  normplot(tdata);

  figure;
  hist(tdata, 100);
end

function [ xp, lambda, c, Lmax ] = boxcox(x, lambda)
  i = find(isfinite(x));
  x = x(i);
  n = length(x);

  xmin = min(x);

  if xmin <= 0
    c = abs(xmin) + 1;
    x = x + c;
  else
    c = 0;
  end

  xmean = mean(x);
  xstd = std(x);

  if nargin < 2
    lambda = fminbnd(@boxcoxf, -10, 10, [], x);
    % lambda = fminsearch(@(l) boxcoxf(l, x), 0);
  end

  if abs(lambda) > eps
    xp = (x.^lambda - 1) / lambda;
  else
    xp = log(x);
  end

  Lmax = -((n - 1) / 2) * log(var(xp)) + ...
    (lambda - 1) * ((n - 1) / n) * sum(log(x));
end

function L = boxcoxf(lambda, x)
  n = length(x);

  if abs(lambda) > eps
    xp = (x.^lambda - 1) / lambda;
  else
    xp = log(x);
  end

  L = -((n - 1) / 2) * log(var(xp)) + ...
    (lambda - 1) * ((n - 1) / n) * sum(log(x));

  L = -L;
end
