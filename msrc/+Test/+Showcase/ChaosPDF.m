init;

c = Config('steps', 100);
c.display();

step = round(c.steps / 3);

ch = Test.constructChaos(c);
[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

for i = 1:c.cores
  Utils.plotHistogram(cRaw(:, i, step), 'unbounded');
  title(sprintf('Chaos: Probability Density of Core %d', i));
end
