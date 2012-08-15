init;

c = Config('steps', 100);
display(c);

orderSet = [ 1 2 3 4 5 ];
sampleSet = [ 10^2 10^3 10^4 10^5 ];

pick = [ 4 4 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

errorExp = zeros(orderCount, sampleCount);
errorVar = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo(c);

fprintf('%5s%15s%15s', 'Order', 'NRMSE(ExpEm)', 'NRMSE(VarEm)');
for k = 1:2
  for i = 1:sampleCount
    fprintf('%15s', sprintf('MC %.1e', sampleSet(i)));
  end
end
fprintf('\n');

mEXP = {};
mVAR = {};

for i = 1:length(orderSet)
  c.tune('constructionMethod.chaosOrder', orderSet(i));

  fprintf('%5d', c.constructionMethod.chaosOrder);

  %% Temperature analysis with Polynomial Chaos.
  %

  ch = Test.constructChaos(c);
  [ cexp, cvar, ~, cexpEm, cvarEm ] = Test.sampleChaos(ch, c);

  fprintf('%15.2f%15.2f', ...
    Stats.NRMSE(cexp, cexpEm) * 100, ...
    Stats.NRMSE(cvar, cvarEm) * 100);

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

    errorExp(i, j) = Stats.NRMSE(mexp, cexp) * 100;
    errorVar(i, j) = Stats.NRMSE(mvar, cvar) * 100;

    if i == pick(1) && j == pick(2)
      mExp = mexp;
      mVar = mvar;
      cExp = cexp;
      cVar = cvar;
    end
  end

  for j = 1:length(sampleSet)
    fprintf('%15.2f', errorExp(i, j));
  end

  for j = 1:length(sampleSet)
    fprintf('%15.2f', errorVar(i, j));
  end

  fprintf('\n');
end

if all(pick == 0), return; end

order = orderSet(pick(1));
samples = sampleSet(pick(2));

time = c.timeLine;

mf = Stats.drawEvolution(time, mExp, mVar);
title(sprintf('%d-sample Monte Carlo', samples));

cf = Stats.drawEvolution(time, cExp, cVar);
title(sprintf('%d-order Polynomial Chaos', order));

Utils.evenScale(mf, cf);

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Error of Expectation');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, (mExp(i, :) - cExp(i, :)) ./ mExp(i, :) * 100, 'Color', color);
end
xlabel('Time, s');
ylabel('Error, %');

subplot(2, 1, 2);
title('Error of Variance');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, (mVar(i, :) - cVar(i, :)) ./ mVar(i, :) * 100, 'Color', color);
end
xlabel('Time, s');
ylabel('Error, %');

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

ExpNRMSE = Stats.NRMSE(mExp, cExp) * 100;
VarNRMSE = Stats.NRMSE(mVar, cVar) * 100;

fprintf('NRMSE(Exp): %.2f%%\n', ExpNRMSE);
fprintf('NRMSE(Var): %.2f%%\n', VarNRMSE);
