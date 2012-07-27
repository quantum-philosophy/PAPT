init;

c = Test.config('cores', 2);
c.display();

hs = HotSpot.Piecewise(c.hotspotSet{:});

t = tic;
T = hs.solve(c.dynamicPower, zeros(hs.sdim, 1));
time = toc(t);
fprintf('Simulation time with Piecewise: %.2f s\n', time);

T = Utils.toCelsius(T);

timeLine = c.timeLine;

figure;
title(sprintf('Solution with Piecewise (%.2f s)', time));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(timeLine, T(i, :), 'Color', color);
end
