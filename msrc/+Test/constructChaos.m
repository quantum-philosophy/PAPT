function ch = constructChaos(c)
  %
  % Construct an appropriate solver.
  %
  ch = HotSpot.Chaos(c.hotspotArguments{:}, c.chaosArguments{:});
  ch.configureLeakage(c.dynamicPower);
end
