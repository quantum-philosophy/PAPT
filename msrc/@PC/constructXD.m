function psiXD = constructXD(x, countXD)
  dimension = length(x);

  %
  % Determine the order of a 1D polynomial.
  %
  check = 0;
  count1D = 0;
  while (check < countXD)
    count1D = count1D + 1;
    check = check + nmultichoosek(dimension, count1D);
  end

  %
  % Create a 1D polynomial.
  %
  psi1D(1, :) = PC.construct1D(x(1), count1D);

  %
  % If there is only one stochastic dimension,
  % we do not need to do anything else.
  %
  if dimension == 1
    psiXD = psi1D;
  else
    %
    % Clone the first 1D polynomial.
    %
    for i = 2:dimension
      psi1D(i, :) = subs(psi1D(1, :), x(1), x(i));
    end

    limit = ones(1, dimension) * (count1D - 1);
    index = zeros(1, dimension);

    done = zeros(1, countXD);

    while (any(limit - index))
      newindex = ndind(index) + 1;

      if newindex <= countXD
        psiXD(newindex) = psi1D(1, index(1) + 1);
        for i = 2:dimension
          psiXD(newindex) = psiXD(newindex) * psi1D(i, index(i) + 1);
        end
        psiXD(newindex) = simplify(expand(psiXD(newindex)));

        done(newindex) = 1;
        if all(done), break; end
      end

      index = genincr(index, limit);
    end
  end
end
