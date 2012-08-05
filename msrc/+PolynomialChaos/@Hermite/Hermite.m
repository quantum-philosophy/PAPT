classdef Hermite < PolynomialChaos.Base
  methods
    function pc = Hermite(dims, order, method)
      pc = pc@PolynomialChaos.Base(dims, order, method);
    end
  end

  methods (Access = 'protected')
    function method = prepare(pc, method)
      if nargin < 2, method = struct(); end

      method = Utils.merge(struct(...
        'quadratureName', 'GaussHermiteProbabilists', ...
        'quadratureType', 'Adaptive'), method);

      method.chaosName = 'Hermite';
    end

    function psi = construct1D(pc, x, order)
      psi(1) = sympoly(1);

      for i = 2:(order + 1)
        psi(i) = x * psi(i - 1) - diff(psi(i - 1), x);
      end
    end
  end
end
