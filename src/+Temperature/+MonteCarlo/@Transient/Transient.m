classdef Transient < Temperature.Numerical.Transient
  properties (SetAccess = 'protected')
    process
  end

  methods
    function this = Transient(varargin)
      options = Options(varargin{:});

      this = this@Temperature.Numerical.Transient(options);
      this.process = ProcessVariation.Normal(options, 'threshold', 1);
    end

    function [ Texp, output ] = compute(this, Pdyn, varargin)
      options = Options(varargin{:});
      verbose = options.get('verbose', false);

      sampleCount = options.get('sampleCount', 1e3);

      filename = options.get('filename', []);
      if isempty(filename)
        filename = sprintf('MonteCarlo_%s.mat', ...
        DataHash({ Pdyn, Utils.toString(this.leakage), ...
          Utils.toString(this.process), sampleCount }));
      end

      if File.exist(filename)
        if verbose
          fprintf('Monte Carlo: using cached data in "%s"...\n', filename);
        end
        load(filename);
      else
        if verbose
          fprintf('Monte Carlo: running %d simulations...\n', sampleCount);
        end

        [ processorCount, stepCount ] = size(Pdyn);

        L = transpose(this.process.sample(sampleCount));

        Tdata = zeros(processorCount, stepCount, sampleCount);

        tic;
        for i = 1:sampleCount
          Tdata(:, :, i) = this.computeWithLeakage( ...
            Pdyn, varargin{:}, 'L', L(:, i));
        end
        time = toc;

        Texp = mean(Tdata, 3);
        Tvar = var(Tdata, [], 3);

        save(filename, 'Texp', 'Tvar', 'Tdata', 'time', '-v7.3');
      end

      if verbose
        fprintf('Monte Carlo: simulation time %.2f s (%d samples).\n', ...
          time, sampleCount);
      end

      Tdata = permute(Tdata, [ 3 1 2 ]);

      output.Tvar = Tvar;
      output.Tdata = Tdata;
      output.time = time;
    end

    function Tdata = evaluate(this, Pdyn, rvs, varargin)
      [ processorCount, stepCount ] = size(Pdyn);

      sampleCount = size(rvs, 1);
      L = transpose(this.process.evaluate(rvs));

      Tdata = zeros(processorCount, stepCount, sampleCount);

      for i = 1:sampleCount
        Tdata(:, :, i) = this.computeWithLeakage( ...
          Pdyn, varargin{:}, 'L', L(:, i));
      end

      Tdata = permute(Tdata, [ 3 1 2 ]);
    end
  end
end
