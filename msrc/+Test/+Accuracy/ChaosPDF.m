init;

c = Config('steps', 100);
display(c);

orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2 10^3 10^4 10^5 ];

pick = [ 0 0 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

error = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo(c);

fprintf('%15s', 'Order');
for i = 1:sampleCount
  fprintf('%15s', sprintf('MC %.1e', sampleSet(i)));
end
fprintf('\n');

mRAW = {};

for i = 1:length(orderSet)
  c.tune('polynomialOrder', orderSet(i));

  fprintf('%15d', c.polynomialOrder);

  %% Temperature analysis with Polynomial Chaos.
  %

  ch = Test.constructChaos(c);
  [ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

  for j = 1:length(sampleSet);
    c.tune('monteCarloSamples', sampleSet(j));

    %% Temperature analysis with Monte Carlo.
    %

    if i == 1
      [ ~, ~, mRAW{j} ] = Test.sampleMonteCarlo(mc, c);
    end

    mRaw = mRAW{j};

    if orderSet(i) == pick(1) && sampleSet(j) == pick(2)
      error(i, j) = Utils.comparePDF(mRaw, cRaw, ...
        c.timeLine, { 'Monte Carlo', 'Polynomial Chaos' }) * 100;
    else
      error(i, j) = Utils.comparePDF(mRaw, cRaw) * 100;
    end

    fprintf('%15.2f', error(i, j));
  end

  fprintf('\n');
end
