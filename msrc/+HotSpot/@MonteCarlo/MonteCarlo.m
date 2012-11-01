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

    function [ Texp, Tvar, Tdata ] = computeWithLeakageInParallel( ...
      this, Pdyn, leakage)

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

      Tdata = zeros(processorCount, stepCount, sampleCount);

      parfor i = 1:sampleCount
        Tdata(:, :, i) = this.computeWithLeakage( ...
          Pdyn, leakage, Lnom + rvMap * rvs(:, i));
      end

      Texp = mean(Tdata, 3);
      Tvar = var(Tdata, [], 3);

      save(filename, 'Texp', 'Tvar', 'Tdata', '-v7.3');
    end
  end
end
