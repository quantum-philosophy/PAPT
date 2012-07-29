function [ Exp, Var, Raw, Time ] = sampleChaos(c)
  hs = HotSpot.Chaos(c.floorplan, c.hotspotConfig, c.hotspotLine, c.order);

  t = tic;

  [ Exp, Var, Trace ] = hs.solve(c.dynamicPower);
  Raw = hs.pc.evaluateSet(Trace, normrnd(0, 1, hs.pc.sdim, c.samples));

  Time = toc(t);

  Exp = Utils.toCelsius(Exp);
  Raw = Utils.toCelsius(Raw);
end
