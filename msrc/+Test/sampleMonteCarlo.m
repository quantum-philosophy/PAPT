function [ exp, var, raw ] = sampleMonteCarlo(one, two)
  if isa(one, 'char')
    mc = Test.constructMonteCarlo(one, two);
    samples = two.samples;
  else
    mc = one;
    samples = two;
  end

  [ exp, var, raw ] = mc.sampleWithReplacement(samples);
  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
