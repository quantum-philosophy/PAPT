function plotDensity(filename)
  setup;
  close all;

  alpha = 0.05;
  mctol = 1e-3;
  shrunkStepCount = 10;
  maximalSampleCount = 1e4;

  load(filename);
  DATA = Tdata;

  [ processorCount, stepCount, sampleCount ] = size(DATA);
  fprintf('Processor count: %d\n', processorCount);
  fprintf('Step count:      %d (%.2f s)\n', stepCount, stepCount * 1e-3);
  fprintf('Sample count:    %d\n', sampleCount);

  sampleCount = min(sampleCount, maximalSampleCount);
  DATA = Utils.toCelsius(DATA(:, :, 1:sampleCount));

  startStep = floor(1 + stepCount / 2 - shrunkStepCount / 2);
  stepCount = shrunkStepCount;
  endStep = startStep + stepCount - 1;
  DATA = DATA(:, startStep:endStep, :);

  k = floor(stepCount / 2);

  fprintf('Considered sample count: %d\n', sampleCount);
  fprintf('Considered step range:   [ %d, %d ]\n', startStep, endStep);
  fprintf('Considered time moment:  %.2f s\n', (startStep + k) * 1e-3);

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

  fprintf('%5s %10s %10s %10s %5s %10s %10s %10s\n', ...
    'H1', 'p1', 'jbstat', 'critval', 'H2', 'p2', 'jbstat', 'critval');
  for i = 1:stepCount
    for j = 1:processorCount
      data = squeeze(DATA(j, i, :));

      if i == k
        [ pdf, x ] = ksdensity(data, 'npoints', sampleCount);
        line(x, pdf, 'Color', Color.pick(j), 'LineWidth', 1);
        labels{end + 1} = sprintf('PE %d', j);
      end

      % continue;

      tdata = boxcox(data);

      [ H1, p1, jbstat1, critval1 ] = jbtest(data, alpha, mctol);
      [ H2, p2, jbstat2, critval2 ] = jbtest(tdata, alpha, mctol);

      fprintf('%5d %10.2e %10.2e %10.2e %5d %10.2e %10.2e %10.2e\n', ...
        H1, p1, jbstat1, critval1, H2, p2, jbstat2, critval2);

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

  for i = 1:processorCount
    figure;

    data = squeeze(DATA(i, k, :));
    subplot(1, 3, 1);
    normplot(data);
    Plot.title('Original');
    Plot.label('Temperature', 'Probability');

    tdata = log(data);
    subplot(1, 3, 2);
    normplot(tdata);
    Plot.title('Log');
    Plot.label('ln(Temperature)', 'Probability');

    tdata = boxcox(data);
    subplot(1, 3, 3);
    normplot(tdata);
    Plot.title('Box-Cox');
    Plot.label('BoxCox(Temperature)', 'Probability');

    figure;
    hist(tdata, 100);
    Plot.title('Box-Cox');
  end
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
