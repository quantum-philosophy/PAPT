init;

c = Config();
c.display();

hs = HotSpot.Analytic(c.hotspotArguments{:});
hs.configureLeakage(c.dynamicPower);

t = tic;
[ T, leakagePower ] = hs.solve(c.dynamicPower, zeros(hs.sdim, 1));
time = toc(t);
fprintf('Simulation time with Analytic: %.2f s\n', time);

T = Utils.toCelsius(T);

timeLine = c.timeLine;

figure;

subplot(2, 1, 1);
title(sprintf('Solution with Analytic (%.2f s)', time));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(timeLine, T(i, :), 'Color', color);
end
xlabel('Time, s');
ylabel('Temperature, C');

subplot(2, 1, 2);
title('Dynamic and Leakage Power');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(timeLine, c.dynamicPower(i, :), 'Color', color);
  line(timeLine, leakagePower(i, :), 'Color', color, 'LineStyle', '--');
end
xlabel('Time, s');
ylabel('Power, W');
