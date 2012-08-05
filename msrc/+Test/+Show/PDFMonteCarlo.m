init;

c = Config();
c.display();

step = round(c.steps / 2);

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Utils.plotHistogram(mRaw(:, i, step), 'unbounded');
  title(sprintf('Monte Carlo: Probability Density of Core %d', i));
end
