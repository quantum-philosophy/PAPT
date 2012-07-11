function psi = construct1D(x, terms)
  psi(1) = ipoly(1);

  for i = 2:terms
    psi(i) = x * psi(i - 1) - diff(psi(i - 1), x);
  end
end
