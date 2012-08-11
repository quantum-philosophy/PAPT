function evenScale(h1, h2)
  lim1 = get(h1, 'ylim');
  lim2 = get(h2, 'ylim');

  lim0 = [ min(lim1(1), lim2(1)) max(lim1(2), lim2(2)) ];

  set(h1, 'ylim', lim0);
  set(h2, 'ylim', lim0);
end
