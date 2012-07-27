init;

c = Test.config('cores', 2);
c.display();

hs = HotSpot.Chaos(c.hotspotSet{:}, c.order);

t = tic;
[ exp, var ] = hs.solve(c.dynamicPower);
time = toc(t);
fprintf('Simulation time with Chaos: %.2f s\n', time);

Utils.plotExpStd(c.timeLine, Utils.toCelsius(exp), var);
title(sprintf('Solution with Polynomial Chaos (%.2f s)', time));
