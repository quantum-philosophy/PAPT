function psi = construct1D(x, order)
  psi(1) = sym(1);

  for i = 2:(order + 1)
    psi(i) = x * psi(i - 1) - diff(psi(i - 1), x);
  end
end
