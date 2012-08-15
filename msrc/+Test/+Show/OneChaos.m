init;

c = Config();
c.display();

hs = HotSpot.Chaos(c.hotspotArguments{:}, c.chaosArguments{:});
hs.configureLeakage(c.dynamicPower);

hs.display();

t = tic;
[ exp, var ] = hs.solve(c.dynamicPower);
time = toc(t);
fprintf('Simulation time with Chaos: %.2f s\n', time);

Stats.drawEvolution(c.timeLine, Utils.toCelsius(exp), var);
title(sprintf('Solution with Polynomial Chaos (%.2f s)', time));
