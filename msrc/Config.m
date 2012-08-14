classdef Config < handle
  %
  % Monte Carlo.
  %
  properties (SetAccess = 'private')
    %
    % The MC method used to assess the PC expansion.
    %
    assessmentMethod = 'Kutta';

    %
    % The total number of loaded samples of the MC method.
    %
    monteCarloTotalSamples = 10^5;

    %
    % The number of samples to use when sampling a PC expansion.
    %
    monteCarloSamples = 10^4;
  end

  %
  % Polynomial chaos expansion.
  %
  properties (SetAccess = 'private')
    %
    % Polynomial chaos and quadrature methods.
    %
    constructionMethod = struct( ...
      'chaosName', 'Hermite', ...
      'chaosType', 'TotalOrder', ...
      'chaosOrder', 4, ...
      ...
      'quadratureName', 'GaussHermiteProbabilists', ...
      'quadratureType', 'Adaptive', ...
      'quadratureOrder', [], ...
      'quadratureLevel', [], ...
      ...
      'jacobiAlpha', 8.1266, ...
      'jacobiBeta', 8.1266, ...
      'jacobiA', -4, ...
      'jacobiB', 4);

    %
    % The number of samples to use when sampling a PC expansion.
    %
    polynomialChaosSamples = 10^5;
  end

  %
  % General system parameters.
  %
  properties (SetAccess = 'private')
    %
    % The number of cores in the system.
    %
    cores = 4;

    %
    % The time step of system profiles.
    %
    dt = 1e-3;

    %
    % The scaling coefficient of the dynamic power.
    %
    powerScale = 2;

    %
    % The number of steps in system profiles.
    %
    steps
  end

  %
  % Automatically determined.
  %
  properties (SetAccess = 'private')
    %
    % The input parameters of HotSpot.
    %
    floorplan
    powerTrace
    hotspotConfig
    hotspotLine

    %
    % The loaded dynamic power profile.
    %
    dynamicPower
  end

  methods
    function c = Config(varargin)
      for i = 1:floor(length(varargin) / 2)
        name = varargin{(i - 1) * 2 + 1};
        value = varargin{(i - 1) * 2 + 2};
        eval(sprintf('c.%s = value;', name));
      end

      %
      % Determine the input files for HotSpot.
      %
      [ c.floorplan, c.powerTrace, c.hotspotConfig ] = ...
        Utils.resolveTest(c.cores);

      c.hotspotLine = sprintf('sampling_intvl %.2e', c.dt);

      %
      % Load the power profile and adjust it if it is needed.
      %
      c.updateDynamicPower();
    end

    function tune(c, varargin)
      for i = 1:floor(length(varargin) / 2)
        name = varargin{(i - 1) * 2 + 1};
        value = varargin{(i - 1) * 2 + 2};

        switch name
        case 'steps'
          c.steps = value;
          c.updateDynamicPower();
        otherwise
          eval(sprintf('c.%s = value;', name));
        end
      end
    end

    function display(c)
      fprintf('Configuration:\n');
      fprintf('  Processing elements: %d\n', c.cores);
      fprintf('  Profile steps: %d\n', c.steps);
      fprintf('  Time step: %.2e s\n', c.dt);
      fprintf('  Simulated time: %.2f s\n', c.dt * c.steps);
      fprintf('  Assessment method: %s\n', c.assessmentMethod);
      fprintf('  Number of MC samples: %d\n', c.monteCarloSamples);
      fprintf('  Polynomial chaos: %s\n', ...
        PolynomialChaos.methodName(c.constructionMethod));
      fprintf('  Integration method: %s\n', ...
        Quadrature.methodName(c.constructionMethod));
    end

    function args = hotspotArguments(c)
      args = { c.floorplan, c.hotspotConfig, c.hotspotLine };
    end

    function args = chaosArguments(c)
      args = { c.constructionMethod };
    end

    function s = stamp(c, varargin)
      p = struct(...
        'assessmentMethod', c.assessmentMethod, ...
        'cores', c.cores, ...
        'steps', c.steps, ...
        'dt', c.dt);

      for i = 1:floor(length(varargin) / 2)
        name = varargin{(i - 1) * 2 + 1};
        value = varargin{(i - 1) * 2 + 2};
        p.(name) = value;
      end

      s = sprintf('%s_nc%d_ns%d_f%d', p.assessmentMethod, p.cores, p.steps, 1 / p.dt);
    end

    function line = timeLine(c)
      line = ((1:c.steps) - 1) * c.dt;
    end
  end

  methods (Access = 'private')
    function updateDynamicPower(c)
      dynamicPower = c.powerScale * dlmread(c.powerTrace, '', 1, 0)';

      if ~isempty(c.steps)
        c.dynamicPower = Utils.replicate(dynamicPower, c.steps);
      else
        c.dynamicPower = dynamicPower;
        c.steps = size(c.dynamicPower, 2);
      end

      assert(size(c.dynamicPower, 1) == c.cores, ...
        'The number of cores is invalid.');
    end
  end
end
