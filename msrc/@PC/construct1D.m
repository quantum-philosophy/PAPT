function psi = construct1D(x, count)
  psi(1) = sym(1);

  for i = 2:count
    psi(i) = x * psi(i - 1) - diff(psi(i - 1), x);
  end
end
