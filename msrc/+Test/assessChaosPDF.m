init;

c = Test.config('steps', 100, 'samples', 10^2);
display(c);

rounds = 10;
orderSet = [ 1 2 3 4 5 6 7 8 9 10 ];
sampleSet = [ 10^2, 10^3, 10^4 ];

orderCount = length(orderSet);
sampleCount = length(sampleSet);

error = zeros(orderCount, sampleCount);

mc = Test.constructMonteCarlo('Kutta', c, max(sampleSet));
ch = Test.constructChaos(c);

fprintf('%15s', 'Order');
for i = 1:sampleCount
  fprintf('%15s', sprintf('MC %.1e', sampleSet(i)));
end
fprintf('\n');

for i = 1:length(orderSet)
  c.order = orderSet(i);

  fprintf('%15d', c.order);

  %% Temperature analysis with Polynomial Chaos.
  %

  [ ~, ~, cRaw ] = Test.sampleChaos(ch, max(sampleSet));

  for j = 1:length(sampleSet);
    c.samples = sampleSet(j);

    %% Temperature analysis with Monte Carlo.
    %

    [ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, sampleSet(j));

    % core = 1;
    % step = round(c.steps / 2);
    % Utils.compareSmooth(mRaw(:, core, step), ...
    %   cRaw(:, core, step), { 'Kutta', 'Chaos' });

    error(i, j) = mean(mean(Utils.comparePDF(mRaw, cRaw), 2)) * 100;

    fprintf('%15.2f', error(i, j));
  end

  fprintf('\n');
end
