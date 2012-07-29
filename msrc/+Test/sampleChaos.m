function [ exp, var, raw ] = sampleChaos(c, samples)
  hs = HotSpot.Chaos(c.hotspotSet{:}, c.order);

  [ exp, var, trace ] = hs.solve(c.dynamicPower);
  exp = Utils.toCelsius(exp);

  if nargout < 3, return; end

  raw = hs.pc.evaluateSet(trace, normrnd(0, 1, hs.pc.sdim, samples));
  raw = Utils.toCelsius(raw);
end
