function [ exp, var, raw ] = sampleChaos(ch, c)
  [ exp, var, trace ] = ch.solve(c.dynamicPower);
  exp = Utils.toCelsius(exp);

  if nargout < 3, return; end

  raw = ch.pc.evaluateSet(trace, normrnd(0, 1, ch.pc.sdim, c.polynomialChaosSamples));
  raw = Utils.toCelsius(raw);
end
