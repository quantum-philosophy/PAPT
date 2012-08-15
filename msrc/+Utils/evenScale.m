function evenScale(one, two)
  one = gca(one);
  two = gca(two);

  oneLimit = get(one, 'ylim');
  twoLimit = get(two, 'ylim');

  zeroLimit = [ min(oneLimit(1), twoLimit(1)) max(oneLimit(2), twoLimit(2)) ];

  set(one, 'ylim', zeroLimit);
  set(two, 'ylim', zeroLimit);
end
