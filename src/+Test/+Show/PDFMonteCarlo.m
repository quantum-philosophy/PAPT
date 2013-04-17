init;

c = Config('steps', 100);
c.display();

step = round(c.steps / 2);

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Stats.observe(mRaw(:, i, step), 'method', 'histogram', 'range', 'unbounded');
  title([ 'Probability Density of Monte Carlo, Core ', num2str(i) ]);
end
