function plotDensity
  setup;

  load('MC.mat');
  DATA = Utils.toCelsius(Tdata2);

  [ sampleCount, processorCount, stepCount ] = size(DATA);

  shrunkStepCount = stepCount;

  startStep = 1 + stepCount / 2 - shrunkStepCount / 2;
  stepCount = shrunkStepCount;

  k = floor(stepCount / 2);

  DATA = DATA(:, :, startStep:(startStep + stepCount - 1));

  figure;

  x = (startStep:(startStep + stepCount - 1))' * 1e-3;

  E = squeeze(mean(DATA, 1));
  V = squeeze(var(DATA, [], 1));

  index = 1:4;

  labels = {};

  for i = index
    I = find(index == i);
    line(x, E(i, :), ...
      'Color', Color.pick(I), 'LineWidth', 1.5);
    line(x, E(i, :) + sqrt(V(i, :)), ...
      'Color', Color.pick(I), 'LineStyle', '--');
    labels{end + 1} = sprintf('PE %d: expectation', I);
    labels{end + 1} = sprintf('PE %d: deviation', I);
  end

  Plot.label('Time, s', 'Temperature, C');
  Plot.legend(labels{:});

  xlim([ x(1), x(end) ]);

  figure;

  pointCount = sampleCount / 10;

  passed1 = 0;
  passed2 = 0;

  labels = {};

  for i = 1:stepCount
    for j = index
      I = find(index == j);

      data = DATA(:, j, i);

      if i == k
        [ pdf, x ] = ksdensity(data, 'npoints', pointCount);
        line(x, pdf, 'Color', Color.pick(I), 'LineWidth', 1);
        labels{end + 1} = sprintf('PE %d', I);
      end

      continue;

      [ tdata, lambda ] = boxcox(data);

      H1 = jbtest(data, 0.001);
      H2 = jbtest(tdata, 0.001);

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

  figure;

  data = DATA(:, index(1), k);
  normplot(data);

  figure;

  [ tdata, lambda ] = boxcox(data);

  m = mean(tdata);
  s = std(tdata);

  tdata = (tdata - m) / s;

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
