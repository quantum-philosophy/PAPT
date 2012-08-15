init;

c = Config('steps', 100);
c.display();

step = round(c.steps / 2);

ch = Test.constructChaos(c);
[ ~, ~, cRaw ] = Test.sampleChaos(ch, c);

mc = Test.constructMonteCarlo(c);
[ ~, ~, mRaw ] = Test.sampleMonteCarlo(mc, c);

for i = 1:c.cores
  Stats.compare(mRaw(:, i, step),  cRaw(:, i, step), ...
    'method', 'histogram', 'draw', true);
  title([ 'Probability Density of Core ', num2str(i) ]);
  legend(c.assessmentMethod, 'Chaos');
end
