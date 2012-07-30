function [ exp, var, raw ] = bootstrapMonteCarlo(one, two, rounds)
  if isa(one, 'char')
    mc = Test.constructMonteCarlo(one, two);
    samples = two.samples;
  else
    mc = one;
    samples = two;
  end

  if nargin < 3, rounds = samples; end

  h = ibar('Bootstrapping Monte Carlo: round %d out of %d.', rounds);

  exp = 0;
  var = 0;
  raw = 0;

  I = [];

  for i = 1:rounds
    [ e, v, r, I ] = mc.sampleAccumulative(samples, I);

    exp = exp + e;
    var = var + v;
    raw = raw + r;

    increase(h);
  end

  exp = exp / rounds;
  var = var / rounds;
  raw = raw / rounds;

  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
