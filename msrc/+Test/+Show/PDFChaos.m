init;

c = Config('steps', 100);
c.display();

step = round(c.steps / 2);

ch = Test.constructChaos(c);
[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

for i = 1:c.cores
  Stats.observe(cRaw(:, i, step), 'method', 'histogram', 'range', 'unbounded');
  title([ 'Probability Density of Chaos, Core ', num2str(i) ]);
end
