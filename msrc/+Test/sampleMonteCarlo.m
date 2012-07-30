function [ exp, var, raw ] = sampleMonteCarlo(mc, c)
  [ exp, var, raw ] = mc.sampleSequential(c.monteCarloSamples);
  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
