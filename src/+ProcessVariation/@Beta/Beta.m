classdef Beta < ProcessVariation.Base
  methods
    function this = Beta(varargin)
      this = this@ProcessVariation.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function variables = distribute(this, variance, correlation, options)
      dimensionCount = length(variance);

      distributions = cell(dimensionCount, 1);

      for i = 1:dimensionCount
        sigma = sqrt(variance(i));

        a = -4 * sigma;
        b =  4 * sigma;

        if i > 1 && variance(i) == variance(i - 1)
          %
          % Use the same parameter.
          %
        else
          param = Utils.fitBetaToNormal('sigma', sigma, ...
            'fitRange', [ a, b ], 'paramRange', [ 1, 20 ]);
        end

        distributions{i} = ProbabilityDistribution.Beta( ...
          'alpha', param, 'beta', param, 'a', a, 'b', b);
      end

      variables = RandomVariables.Heterogeneous( ...
        'distributions', distributions, 'correlation', correlation);
    end
  end
end
