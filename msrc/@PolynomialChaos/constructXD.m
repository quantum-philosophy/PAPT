function psiXD = constructXD(x, termsXD)
  sdim = length(x);

  %
  % Determine the order of a 1D polynomial.
  %
  check = 0;
  terms1D = 0;
  while (check < termsXD)
    terms1D = terms1D + 1;
    check = check + nmultichoosek(sdim, terms1D);
  end

  %
  % Create a 1D polynomial.
  %
  psi1D(1, :) = PolynomialChaos.construct1D(x(1), terms1D);

  %
  % If there is only one stochastic dimension,
  % we do not need to do anything else.
  %
  if sdim == 1
    psiXD = psi1D;
  else
    %
    % Clone the first 1D polynomial.
    %
    for i = 2:sdim
      psi1D(i, :) = subs(psi1D(1, :), x(1), x(i));
    end

    limit = ones(1, sdim) * (terms1D - 1);
    index = zeros(1, sdim);

    done = zeros(1, termsXD);

    found = 0;
    h = waitbar(found / termsXD, ...
      sprintf('Polynomial Chaos, terms %d/%d.', found, termsXD));

    while (any(limit - index))
      newindex = ndind(index) + 1;

      if newindex <= termsXD
        psiXD(newindex) = psi1D(1, index(1) + 1);
        for i = 2:sdim
          psiXD(newindex) = psiXD(newindex) * psi1D(i, index(i) + 1);
        end

        done(newindex) = 1;

        found = found + 1;
        waitbar(found / termsXD, h, ...
          sprintf('Polynomial Chaos, terms %d/%d.', found, termsXD));

        if all(done), break; end
      end

      index = genincr(index, limit);
    end

    close(h);
  end
end
