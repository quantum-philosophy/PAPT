init;

c = Test.config();
c.adjustPowerSteps(100);
c.display();

divide = 2;

X = [ 1 2 3 4 5 6 ];

pick = length(X);

%% Temperature analysis with Monte Carlo.
%

[ mExp, mVar, ~, mTime ] = Test.computeKutta(c);

fprintf('MC simulation time: %.2f s\n', mTime);

%% Polynomial Chaos expansion.
%

Pdyn = zeros(c.cores, c.steps * divide);

I = ((1:c.steps) - 1) * divide + 1;
for i = 1:divide
  Pdyn(:, I + i - 1) = c.dynamicPower;
end
I = I + (divide - 1);

fprintf('Extended number of steps: %d\n', c.steps * divide);
fprintf('Modified sampling interval: %.2e s\n', c.dt / divide);

fprintf('%15s%15s%15s%15s\n', 'Order', 'Time, s', 'NRMSE(Exp), %', 'NRMSE(Var), %');

count = length(X);

errorExp = zeros(count, 1);
errorVar = zeros(count, 1);

c.adjustSamplingInterval(c.dt / divide);

for i = 1:count
  order = X(i);
  hs = HotSpot.Chaos(c.floorplan, c.hotspotConfig, c.hotspotLine, order);
  assert(hs.dt == c.dt, 'The sampling interval is invalid.');

  t = tic;
  [ exp, var ] = hs.solve(Pdyn);
  time = toc(t);

  exp = Utils.toCelsius(exp(:, I));
  var = var(:, I);

  errorExp(i) = Utils.NRMSE(mExp, exp) * 100;
  errorVar(i) = Utils.NRMSE(mVar, var) * 100;

  fprintf('%15d%15.2f%15.2f%15.2f\n', order, time, errorExp(i), errorVar(i));

  if pick == i
    cExp = exp;
    cVar = var;
    cTime = time;
  end
end

time = c.timeLine;

mStd = sqrt(mVar);
cStd = sqrt(cVar);

mf = figure;
title(sprintf('%d-sample Monte Carlo (%.2f s)', c.samples, mTime));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mExp(i, :), 'Color', color);
  line(time, mExp(i, :) + mStd(i, :), 'Color', color, 'LineStyle', '--');
  line(time, mExp(i, :) - mStd(i, :), 'Color', color, 'LineStyle', '--');
end

cf = figure;
title(sprintf('%d-order Polynomial Chaos (%.2f s)', X(pick), cTime));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, cExp(i, :), 'Color', color);
  line(time, cExp(i, :) + cStd(i, :), 'Color', color, 'LineStyle', '--');
  line(time, cExp(i, :) - cStd(i, :), 'Color', color, 'LineStyle', '--');
end

%% Make the plots match each other.
%

lim1 = ylim;
figure(mf);
lim2 = ylim;
lim = [ min(lim1(1), lim2(1)) max(lim1(2), lim2(2)) ];
ylim(lim);
figure(cf);
ylim(lim);

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Difference of Expectations');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mExp(i, :) - cExp(i, :), 'Color', color);
end

subplot(2, 1, 2);
title('Difference of Variances');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mVar(i, :) - cVar(i, :), 'Color', color);
end

figure;
subplot(2, 1, 1);
title('Convergence of Expectations');
color = Utils.pickColor(1);
line(X, errorExp, 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE of Expectation, %');

subplot(2, 1, 2);
title('Convergence of Variances');
color = Utils.pickColor(1);
line(X, errorVar, 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE of Variance, %');
