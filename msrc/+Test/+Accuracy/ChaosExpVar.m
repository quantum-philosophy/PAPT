init;

c = Config('steps', 100);
display(c);

orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2, 10^3 10^4 10^5 ];

pick = [ 4, 10^4 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

errorExp = zeros(orderCount, sampleCount);
errorVar = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo(c);

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
  c.tune('polynomialOrder', orderSet(i));

  fprintf('%15d', c.polynomialOrder);

  %% Temperature analysis with Polynomial Chaos.
  %

  ch = Test.constructChaos(c);
  [ cexp, cvar ] = Test.sampleChaos(ch, c);

  for j = 1:length(sampleSet);
    c.tune('monteCarloSamples', sampleSet(j));

    %% Temperature analysis with Monte Carlo.
    %

    if i == 1
      [ mEXP{j}, mVAR{j} ] = Test.sampleMonteCarlo(mc, c);
    end

    mexp = mEXP{j};
    mvar = mVAR{j};

    %% Comparison of the methods.
    %

    errorExp(i, j) = Utils.NRMSE(mexp, cexp) * 100;
    errorVar(i, j) = Utils.NRMSE(mvar, cvar) * 100;

    if orderSet(i) == pick(1) && sampleSet(j) == pick(2)
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

mf = Utils.plotExpStd(time, mExp, mVar);
title(sprintf('%d-sample Monte Carlo', sampleSet(pick(2))));

cf = Utils.plotExpStd(time, cExp, cVar);
title(sprintf('%d-order Polynomial Chaos', orderSet(pick(1))));

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
title('Difference of Expectation');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mExp(i, :) - cExp(i, :), 'Color', color);
end
xlabel('Time, s');
ylabel('Exp(PC) - Exp(MC), C');

subplot(2, 1, 2);
title('Difference of Variance');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, mVar(i, :) - cVar(i, :), 'Color', color);
end
xlabel('Time, s');
ylabel('Var(PC) - Var(MC), $C^2$');

figure;
subplot(2, 1, 1);
title('Convergence of Expectation');
color = Utils.pickColor(1);
line(orderSet, errorExp(:, pick(2)), 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE(Exp), %');

subplot(2, 1, 2);
title('Convergence of Variance');
color = Utils.pickColor(1);
line(orderSet, errorVar(:, pick(2)), 'Color', color);
xlabel('Polynomial Order');
ylabel('NRMSE(Var), %');