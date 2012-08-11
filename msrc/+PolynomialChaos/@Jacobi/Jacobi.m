classdef Jacobi < PolynomialChaos.Base
  properties (SetAccess = 'private')
    alpha
    beta
    a
    b
  end

  methods
    function pc = Jacobi(varargin)
      pc = pc@PolynomialChaos.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function method = prepare(pc, method)
      method = Utils.merge(struct( ...
        'chaosType', 'TotalOrder', ...
        'quadratureName', 'GaussJacobi', ...
        'quadratureType', 'Adaptive'), method);

      method.chaosName = 'Jacobi';

      pc.alpha = method.jacobiAlpha;
      pc.beta = method.jacobiBeta;
      pc.a = method.jacobiA;
      pc.b = method.jacobiB;
    end

    function psi = construct1D(pc, X, order)
      alpha = pc.alpha;
      beta = pc.beta;

      psi(1) = sympoly(1);
      x = sym(char(X));

      for n = 0:order
        p = ((-1)^n / (2^n * factorial(n))) * ...
          (1 - x)^(-alpha) * (1 + x)^(-beta) * ...
          diff((1 - x)^(n + alpha) * (1 + x)^(n + beta), x, n);

        psi(n + 1) = sympoly.convertsym(expand(simplify(p)));
      end
    end

    function nodes = generateSampleNodes(pc, points)
      nodes = ibetarnd(pc.alpha + 1, pc.beta + 1, -1, 1, pc.sdim, points);
    end
  end
end
