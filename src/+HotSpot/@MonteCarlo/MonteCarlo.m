classdef MonteCarlo < HotSpot.Numeric
  properties (SetAccess = 'protected')
    process

    sampleCount
    filename
    verbose
  end

  methods
    function this = MonteCarlo(varargin)
      options = Options(varargin{:});

      this = this@HotSpot.Numeric(options);

      this.process = ProcessVariation( ...
        'expectation', LeakagePower.Lnom, ...
        'deviation', 0.05 * LeakagePower.Lnom, ...
        'threshold', '1', options);

      this.sampleCount = options.get('sampleCount', 1e3);
      this.filename = options.get('filename', []);
      if options.get('verbose', false)
        this.verbose = @(varargin) fprintf(varargin{:});
      else
        this.verbose = @(varargin) [];
      end
    end

    function [ Texp, Tvar, Tdata ] = computeInParallel(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      sampleCount = this.sampleCount;

      verbose = this.verbose;

      filename = this.filename;
      if isempty(filename)
        filename = sprintf('MonteCarlo_%s.mat', ...
          DataHash({ Pdyn, Utils.toString(leakage), sampleCount }));
      end

      if File.exist(filename)
        verbose('Monte Carlo: using cached data in "%s"...\n', filename);
        load(filename);
        if ~exist('time', 'var'), time = 0; end
      else
        verbose('Monte Carlo: running %d simulations...\n', sampleCount);

        process = this.process;

        rvs = normrnd(0, 1, process.dimensionCount, sampleCount);

        expectation = process.expectation;
        deviation = process.deviation;
        mapping = process.mapping;

        Tdata = zeros(processorCount, stepCount, sampleCount);

        tic;
        parfor i = 1:sampleCount
          Tdata(:, :, i) = this.compute( ...
            Pdyn, leakage, expectation + deviation * mapping * rvs(:, i));
        end
        time = toc;

        Texp = mean(Tdata, 3);
        Tvar = var(Tdata, [], 3);

        save(filename, 'Texp', 'Tvar', 'Tdata', 'time', '-v7.3');
      end
      verbose('Monte Carlo: simulation time %.2f s (%d samples).\n', ...
        time, sampleCount);

      Tdata = permute(Tdata, [ 3 1 2 ]);
    end

   function Tdata = evaluateInParallel(this, Pdyn, leakage, rvs)
      [ processorCount, stepCount ] = size(Pdyn);

      process = this.process;

      rvs = process.expectation + process.deviation * process.mapping * rvs.';
      sampleCount = size(rvs, 2);

      Tdata = zeros(processorCount, stepCount, sampleCount);

      parfor i = 1:sampleCount
        Tdata(:, :, i) = this.compute(Pdyn, leakage, rvs(:, i));
      end

      Tdata = permute(Tdata, [ 3 1 2 ]);
    end
  end
end
