classdef Hermite < PolynomialChaos.Base
  methods
    function pc = Hermite(dims, order, method)
      if nargin < 3
        method = struct( ...
          'quadratureName', 'GaussHermiteProbabilists', ...
          'quadratureType', 'Adaptive');
      end
      method.chaosName = 'Hermite';
      pc = pc@PolynomialChaos.Base(dims, order, method);
    end
  end

  methods (Access = 'protected')
    function psi = construct1D(pc, x, order)
      psi(1) = sympoly(1);

      for i = 2:(order + 1)
        psi(i) = x * psi(i - 1) - diff(psi(i - 1), x);
      end
    end
  end
end
