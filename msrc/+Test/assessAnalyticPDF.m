init;

c = Test.config('steps', 100);
c.display();

[ kExp, kVar, kRaw, kTime ] = Test.sampleKutta(c);
[ aExp, aVar, aRaw, aTime ] = Test.sampleAnalytic(c);

error = Utils.comparePDF(kRaw, aRaw, c.timeLine, { 'Kutta', 'Analytic' });

fprintf('%10s%20s\n', 'Core', 'Average NRMSE, %');
for i = 1:c.cores
  fprintf('%10d%20.2f\n', i, mean(error(i, :)) * 100);
end
