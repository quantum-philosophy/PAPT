init;

c = Config('steps', 100);

step = round(c.steps / 3);

mc = Test.constructMonteCarlo(c);
ch = Test.constructChaos(c);

[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Utils.compareHistogram(mRaw(:, i, step), ...
    cRaw(:, i, step), { c.assessmentMethod, 'Chaos' });
  title(sprintf('Probability Density for Core %d', i));
end
