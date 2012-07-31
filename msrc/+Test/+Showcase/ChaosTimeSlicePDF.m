%% Initialize.
%

init;

c = Config('steps', 100);

step = round(c.steps / 3);

%% Analyse with Polynomial Chaos.
%

ch = Test.constructChaos(c);
[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

%% Analyse with Monte Carlo.
%

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Utils.compareHistogram(mRaw(:, i, step), ...
    cRaw(:, i, step), { c.assessmentMethod, 'Chaos' });
  title(sprintf('Empirical Probability Density of Core %d', i));
end
