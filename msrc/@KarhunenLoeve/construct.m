function [ values, functions ] = construct(this, options)
  tolerance = 1e-10;
  epsilon = 1e-6;

  a = options.domainBoundary;
  c = options.correlationLength;
  d = options.dimension;

  omegas = zeros(1, d);

  even = @(x) c - x * tan(a * x);
  odd  = @(x) x + c * tan(a * x);

  for i = 0:2:(d - 1)
    left  = (pi / (2 * a)) * (i    ) + epsilon;
    right = (pi / (2 * a)) * (i + 1) - epsilon;

    omegas(i + 1) = bisect(left, right, even);

    left  = (pi / (2 * a)) * (i + 1) + epsilon;
    right = (pi / (2 * a)) * (i + 2) - epsilon;

    omegas(i + 2) = bisect(left, right, odd);
  end

  values = zeros(1, d);
  for i = 1:d
    values(i) = (2 * c) / (omegas(i)^2 + c^2);
  end

  functions = cell(1, d);
  for i = 1:d
    omega = omegas(i);

    if mod(i - 1, 2) == 1
      functions{i} = @(x) sin(omega * x) / ...
        sqrt(a - (sin(2 * omega * a)) / (2 * omega));
    else
      functions{i} = @(x) cos(omega * x) / ...
        sqrt(a + (sin(2 * omega * a)) / (2 * omega));
    end
  end
end
