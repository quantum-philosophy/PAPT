init;

c = Config();
display(c);

orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2 10^3 10^4 ];

pick = [ 0 0 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

error = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo(c);

mEXP = cell(sampleCount, 1);
mVAR = cell(sampleCount, 1);
mRAW = cell(sampleCount, 1);

errorExp = zeros(orderCount, sampleCount);
errorVar = zeros(orderCount, sampleCount);
errorPDF = zeros(orderCount, sampleCount);

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
    fprintf('%15s', sprintf('MC %.1e', sampleSet(j)));
  end
end
fprintf('\n');

for i = 1:length(orderSet)
  c.tune('polynomialOrder', orderSet(i));

  fprintf('%5d | ', c.polynomialOrder);

  %% Temperature analysis with Polynomial Chaos.
  %

  ch = Test.constructChaos(c);
  [ cExp, cVar, cRaw ] = Test.sampleChaos(ch, c);

  for j = 1:length(sampleSet);
    c.tune('monteCarloSamples', sampleSet(j));

    %% Temperature analysis with Monte Carlo.
    %

    if i == 1
      [ mEXP{j}, mVAR{j}, mRAW{j} ] = Test.sampleMonteCarlo(mc, c);
    end

    mExp = mEXP{j};
    mVar = mVAR{j};
    mRaw = mRAW{j};

    errorExp(i, j) = Utils.NRMSE(mExp, cExp) * 100;
    errorVar(i, j) = Utils.NRMSE(mVar, cVar) * 100;

    if orderSet(i) == pick(1) && sampleSet(j) == pick(2)
      % core = 1;
      % step = round(c.steps / 2);
      % Utils.compareSmooth(mRaw(:, core, step), ...
      %   cRaw(:, core, step), { 'Kutta', 'Chaos' });

      errorPDF(i, j) = Utils.comparePDF(mRaw, cRaw, ...
        c.timeLine, { 'Monte Carlo', 'Polynomial Chaos' }) * 100;
    else
      errorPDF(i, j) = Utils.comparePDF(mRaw, cRaw) * 100;
    end

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
