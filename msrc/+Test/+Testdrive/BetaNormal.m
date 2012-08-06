clc;

%% The desired normal distribution.
%
mu = 0;
sigma = 1;
times = 4;

max = mu + times * sigma;
min = mu - max;

x = min:0.01:max;

n = normpdf(x, mu, sigma);

%% The approximation with the Beta distribution.
%
if 0
  rvs = mu + sigma * randn(1, 1000);
  I = find(abs(rvs) > (times * sigma));
  rvs(I) = [];
  rvs = (rvs + times * sigma) / (2 * times * sigma);

  [ phat, pci ] = betafit(rvs);

  alpha = phat(1);
else
  f = @(a) Utils.NRMSE(n, ibetapdf(x, a, a, min, max));
  alpha = fminsearch(f, 8);
end

beta = alpha;

fprintf('Alpha, Beta = %.4f\n', alpha);

b = ibetapdf(x, alpha, beta, min, max);

figure;
title('Probability Density Function');

line(x, b, 'Color', Utils.pickColor(1));
line(x, n, 'Color', Utils.pickColor(2));

b = ibetarnd(alpha, beta, min, max, 10^5, 1);
b = ksdensity(b, x);
line(x, b, 'Color', Utils.pickColor(3));

legend('Beta', 'Normal', 'Monte Carlo');

xlabel('x');
ylabel('PDF');
