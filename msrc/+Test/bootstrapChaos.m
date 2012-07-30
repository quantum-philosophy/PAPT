function [ exp, var, raw ] = bootstrapChaos(c, samples, rounds)
  hs = HotSpot.Chaos(c.hotspotSet{:}, c.order);

  h = ibar('Bootstrapping Chaos: round %d out of %d.', rounds);

  raw = 0;

  [ exp, var, trace ] = hs.solve(c.dynamicPower);
  exp = Utils.toCelsius(exp);

  for i = 1:rounds
    r = hs.pc.evaluateSet(trace, normrnd(0, 1, hs.pc.sdim, samples));

    raw = raw + r;

    increase(h);
  end

  raw = raw / rounds;

  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
