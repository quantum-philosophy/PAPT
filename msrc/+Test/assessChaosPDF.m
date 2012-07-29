init;

c = Test.config('steps', 100);
display(c);

X = [ 1 2 3 4 5 6 ];

%% Temperature analysis with Monte Carlo.
%

[ mExp, mVar, mRaw, mTime ] = Test.sampleKutta(c);

%% Temperature analysis with Polynomial Chaos.
%

cores = c.cores;

fprintf('%15s%15s', 'Order', 'Time, s');
labels = {};
for i = 1:cores
  labels{end + 1} = sprintf('NRMSE(%d), %%', i);
end
labels{end + 1} = 'NRMSE, %';
fprintf('%15s', labels{:});
fprintf('\n');

count = length(X);

for i = 1:count
  c.order = X(i);

  [ cExp, cVar, cRaw, cTime ] = Test.sampleChaos(c);

  % core = 1;
  % step = round(c.steps / 2);
  % Utils.compareSmooth(mRaw(:, core, step), ...
  %   cRaw(:, core, step), { 'Kutta', 'Chaos' });

  error = Utils.comparePDF(mRaw, cRaw);
  error = mean(error, 2);

  fprintf('%15d%15.2f', c.order, cTime);
  fprintf('%15.2f', [ error; mean(error) ] * 100);
  fprintf('\n');
end
