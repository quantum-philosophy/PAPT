classdef MonteCarlo < HotSpot.Analytic & ProcessVariation
  properties (SetAccess = 'private')
    sampleCount
  end

  methods
    function this = MonteCarlo(floorplan, config, line, varargin)
      options = Options(varargin{:});

      this = this@HotSpot.Analytic(floorplan, config, line);
      this = this@ProcessVariation(floorplan, options, ...
        'reduction', 'none');

      assert(this.rvCount == this.processorCount + 1);

      this.sampleCount = options.get('sampleCount', 1e3);
    end

    function [ Texp, Tvar ] = computeWithLeakage(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      sampleCount = this.sampleCount;

      filename = sprintf('MonteCarlo_%s.mat', ...
        DataHash({ Pdyn, Utils.toString(leakage), sampleCount }));

      if File.exist(filename)
        load(filename);
        return;
      end

      rvs = normrnd(0, 1, this.rvCount, sampleCount);

      Lnom = leakage.Lnom;
      rvMap = this.rvMap;

      data = zeros(processorCount, stepCount, sampleCount);

      for i = 1:sampleCount
        data(:, :, i) = computeWithLeakage@HotSpot.Analytic( ...
          this, Pdyn, leakage, Lnom + rvMap * rvs(:, i));
      end

      Texp = mean(data, 3);
      Tvar = var(data, [], 3);

      save(filename, 'Texp', 'Tvar', '-v7.3');
    end
  end
end
