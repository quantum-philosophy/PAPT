init;

trials = 10;

c = Config('steps', 100);
mc = Test.constructMonteCarlo(c);

sampleSet = [ 10^2 10^3 10^4 ];
sampleCount = length(sampleSet);

time = c.timeLine;

for i = 1:sampleCount
  samples = sampleSet(i);

  I = 1:samples;

  Exp = zeros(c.cores, c.steps, trials);
  Var = zeros(c.cores, c.steps, trials);

  for j = 1:trials
    [ Exp(:, :, j), Var(:, :, j), ~ ] = mc.sampleWithReplacement(samples, I);
  end

  Exp = var(Exp, 1, 3);
  Var = var(Var, 1, 3);

  figure;

  subplot(2, 1, 1);
  line(time, Exp);
  title(sprintf('Var(Exp) with %d samples', samples));

  subplot(2, 1, 2);
  line(time, Var);
  title(sprintf('Var(Var) with %d samples', samples));
end
