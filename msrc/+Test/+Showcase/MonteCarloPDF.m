init;

c = Config('steps', 100);
c.display();

step = round(c.steps / 3);

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Utils.plotHistogram(mRaw(:, i, step), 'unbounded');
  title(sprintf('Monte Carlo: Probability Density of Core %d', i));
end
