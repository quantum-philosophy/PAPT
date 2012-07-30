init;

c = Test.config('steps', 100, 'samples', 0);
display(c);

rounds = 1;

orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2, 10^3, 10^4, 10^5 ];

pick = [ length(orderSet), length(sampleSet) ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

errorExp = zeros(orderCount, sampleCount);
errorVar = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo('Kutta', c, 1e5);
ch = Test.constructChaos(c);

fprintf('%15s', 'Order');
for k = 1:2
  for i = 1:sampleCount
    fprintf('%15s', sprintf('MC %.1e', sampleSet(i)));
  end
end
fprintf('\n');

mEXP = {};
mVAR = {};

for i = 1:length(orderSet)
  c.order = orderSet(i);

  fprintf('%15d', c.order);

  %% Temperature analysis with Polynomial Chaos.
  %

  [ cexp, cvar ] = Test.sampleChaos(ch);

  for j = 1:length(sampleSet);
    c.samples = sampleSet(j);

    %% Temperature analysis with Monte Carlo.
    %

    if i == 1
      [ mEXP{j}, mVAR{j} ] = Test.bootstrapMonteCarlo(mc, sampleSet(j), rounds);
    end

    mexp = mEXP{j};
    mvar = mVAR{j};

    %% Comparison of the methods.
    %

    errorExp(i, j) = Utils.NRMSE(mexp, cexp) * 100;
    errorVar(i, j) = Utils.NRMSE(mvar, cvar) * 100;

    if i == pick(1) && j == pick(2)
      mExp = mexp;
      mVar = mvar;
      cExp = cexp;
      cVar = cvar;
    end
  end

  for j = 1:length(sampleSet);
    fprintf('%15.2f', errorExp(i, j));
  end

  for j = 1:length(sampleSet);
    fprintf('%15.2f', errorVar(i, j));
  end

  fprintf('\n');
end

time = c.timeLine;

mStd = sqrt(mVar);
cStd = sqrt(cVar);

mf = figure;
title(sprintf('%d-sample Monte Carlo', c.samples));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mExp(i, :), 'Color', color);
  line(time, mExp(i, :) + mStd(i, :), 'Color', color, 'LineStyle', '--');
  line(time, mExp(i, :) - mStd(i, :), 'Color', color, 'LineStyle', '--');
end

cf = figure;
title(sprintf('%d-order Polynomial Chaos', orderSet(pick(1))));
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
line(orderSet, errorExp(:, pick(2)), 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE of Expectation, %');

subplot(2, 1, 2);
title('Convergence of Variances');
color = Utils.pickColor(1);
line(orderSet, errorVar(:, pick(2)), 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE of Variance, %');
