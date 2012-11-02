classdef MonteCarlo < HotSpot.Analytic & ProcessVariation
  properties (SetAccess = 'private')
    sampleCount
    verbose
  end

  methods
    function this = MonteCarlo(floorplan, config, line, varargin)
      options = Options(varargin{:});

      this = this@HotSpot.Analytic(floorplan, config, line);
      this = this@ProcessVariation(floorplan, options, ...
        'reduction', 'none');

      assert(this.rvCount == this.processorCount + 1);

      this.sampleCount = options.get('sampleCount', 1e3);
      if options.get('verbose', false)
        this.verbose = @(varargin) fprintf(varargin{:});
      else
        this.verbose = @(varargin) [];
      end
    end

    function [ Texp, Tvar, Tdata ] = computeWithLeakageInParallel( ...
      this, Pdyn, leakage)

      [ processorCount, stepCount ] = size(Pdyn);
      sampleCount = this.sampleCount;

      verbose = this.verbose;

      filename = sprintf('MonteCarlo_%s.mat', ...
        DataHash({ Pdyn, Utils.toString(leakage), sampleCount }));

      if File.exist(filename)
        verbose('Monte Carlo: using cached data in "%s"...\n', filename);
        load(filename);
        if ~exist('time', 'var'), time = 0; end
      else
        verbose('Monte Carlo: running %d simulations...\n', sampleCount);

        rvs = normrnd(0, 1, this.rvCount, sampleCount);

        Lnom = leakage.Lnom;
        rvMap = this.rvMap;

        Tdata = zeros(processorCount, stepCount, sampleCount);

        tic;
        parfor i = 1:sampleCount
          Tdata(:, :, i) = this.computeWithLeakage( ...
            Pdyn, leakage, Lnom + rvMap * rvs(:, i));
        end
        time = toc;

        Texp = mean(Tdata, 3);
        Tvar = var(Tdata, [], 3);

        save(filename, 'Texp', 'Tvar', 'Tdata', 'time', '-v7.3');
      end

      verbose('Monte Carlo: simulation time %.2f s.\n', time);
    end
  end
end
