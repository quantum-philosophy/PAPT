init;

c = Test.config('steps', 100, 'samples', 0);
display(c);

chaosSamples = 10^5;

orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2, 10^3, 10^4, 10^5 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

error = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo('Kutta', c, 10^5);
ch = Test.constructChaos(c);

mEXP = cell(sampleCount, 1);
mVAR = cell(sampleCount, 1);
mRAW = cell(sampleCount, 1);

errorExp = zeros(orderCount, sampleCount);
errorVar = zeros(orderCount, sampleCount);
errorPDF = zeros(orderCount, sampleCount);

fprintf('\n');

names = { 'NRMSE(PDF)', 'NRMSE(Exp)', 'NRMSE(Var)' };

fprintf('%15s | ', '');
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

fprintf('%15s | ', 'Order');
for i = 1:length(names)
  if i > 1, fprintf(' | '); end

  for j = 1:sampleCount
    fprintf('%15s', sprintf('MC %.1e', sampleSet(j)));
  end
end
fprintf('\n');

for i = 1:length(orderSet)
  c.order = orderSet(i);

  fprintf('%15d | ', c.order);

  %% Temperature analysis with Polynomial Chaos.
  %

  [ cExp, cVar, cRaw ] = Test.sampleChaos(ch, chaosSamples);

  for j = 1:length(sampleSet);
    c.samples = sampleSet(j);

    %% Temperature analysis with Monte Carlo.
    %

    if i == 1
      [ mEXP{j}, mVAR{j}, mRAW{j} ] = Test.sampleMonteCarlo(mc, sampleSet(j));
    end

    mExp = mEXP{j};
    mVar = mVAR{j};
    mRaw = mRAW{j};

    errorExp(i, j) = Utils.NRMSE(mExp, cExp) * 100;
    errorVar(i, j) = Utils.NRMSE(mVar, cVar) * 100;

    % core = 1;
    % step = round(c.steps / 2);
    % Utils.compareSmooth(mRaw(:, core, step), ...
    %   cRaw(:, core, step), { 'Kutta', 'Chaos' });

    errorPDF(i, j) = mean(mean(Utils.comparePDF(mRaw, cRaw), 2)) * 100;

    fprintf('%15.2f', errorPDF(i, j));
  end

  fprintf(' | ');

  for j = 1:length(sampleSet);
    fprintf('%15.2f', errorExp(i, j));
  end

  fprintf(' | ');

  for j = 1:length(sampleSet);
    fprintf('%15.2f', errorVar(i, j));
  end

  fprintf('\n');
end
