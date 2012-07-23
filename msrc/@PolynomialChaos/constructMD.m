function [ psiMD, index ] = constructMD(x, order)
  sdim = length(x);

  index = Utils.computeMultiIndex(ones(1, sdim) * order) + 1;

  %
  % Create a 1D polynomial.
  %
  psi1D(1, :) = PolynomialChaos.construct1D(x(1), order);

  %
  % If there is only one stochastic dimension,
  % we do not need to do anything else.
  %
  if sdim == 1
    psiMD = psi1D(1, :);
    return;
  end

  %
  % Clone the first 1D polynomial.
  %
  for i = 2:sdim
    psi1D(i, :) = subs(psi1D(1, :), x(1), x(i));
  end

  terms = size(index, 1);

  for i = 1:terms
    psiMD(i) = psi1D(1, index(i, 1));
    for j = 2:sdim
      psiMD(i) = psiMD(i) * psi1D(j, index(i, j));
    end
  end
end
