classdef MultiProbabilists < GaussianQuadrature.Probabilists
  methods
    function gq = MultiProbabilists(x, psi, index)
      order = PolynomialChaos.indexToOrder(index);
      gq = gq@GaussianQuadrature.Probabilists(x, psi, order, index);
    end
  end

  methods (Access = 'protected')
    [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(gq, x, psi, order, index)
  end
end
