init;

c = Config();
c.display();

step = round(c.steps / 2);

ch = Test.constructChaos(c);
[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Utils.compareHistogram(mRaw(:, i, step), ...
    cRaw(:, i, step), { c.assessmentMethod, 'Chaos' });
  title(sprintf('Probability Density of Core %d', i));
end
