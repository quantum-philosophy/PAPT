function [ Exp, Var, Raw, Time ] = sample(method, c, rounds)
  if nargin < 2, c = Test.config(); end
  if nargin < 3, rounds = 1; end

  hs = HotSpot.(method)(c.hotspotSet{:});

  Exp = 0;
  Var = 0;
  Raw = 0;
  Time = 0;

  if rounds > 1
    h = ibar([ sprintf('Monte Carlo with %s', method), ': round %d out of %d.' ], rounds);
  end

  for i = 1:rounds
    [ exp, var, raw, time ] = MonteCarlo.perform3D( ...
      @(rvs) hs.solve(c.dynamicPower, rvs), [ hs.sdim, c.cores, c.steps ], ...
      c.samples, c.stamp(method));

    Exp = Exp + exp;
    Var = Var + var;
    Raw = Raw + raw;
    Time = Time + time;

    if rounds > 1, increase(h); end
  end

  Exp = Exp ./ rounds;
  Var = Var ./ rounds;
  Raw = Raw ./ rounds;
  Time = Time ./ rounds;

  Exp = Utils.toCelsius(Exp);
  Raw = Utils.toCelsius(Raw);
end
