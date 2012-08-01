init;

c = Config('cores', 2);
c.display();

hs = HotSpot.Chaos(c.hotspotArguments{:}, c.chaosArguments{:});

t = tic;
[ exp, var ] = hs.solve(c.dynamicPower);
time = toc(t);
fprintf('Simulation time with Chaos: %.2f s\n', time);

Utils.plotExpStd(c.timeLine, Utils.toCelsius(exp), var);
title(sprintf('Solution with Polynomial Chaos (%.2f s)', time));
