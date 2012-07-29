classdef config < handle
  properties
    %
    % HotSpot inputs.
    %
    floorplan
    powerTrace
    hotspotConfig
    hotspotLine

    %
    % Basic parameters.
    %
    cores
    steps
    dt

    %
    % Monte Carlo sampling
    %
    samples

    %
    % Polynomial chaos expansion.
    %
    order

    %
    % Dynamic power profile.
    %
    originalDynamicPower
    dynamicPower
  end

  methods
    function c = config(varargin)
      c.cores = 4;
      c.dt = 1e-3;
      c.samples = 1e4;
      c.order = 4;

      for i = 1:floor(length(varargin) / 2)
        name = lower(varargin{(i - 1) * 2 + 1});
        value = varargin{(i - 1) * 2 + 2};
        c.(name) = value;
      end

      [ c.floorplan, c.powerTrace, c.hotspotConfig ] = ...
        Utils.resolveTest(c.cores);

      c.hotspotLine = sprintf('sampling_intvl %.2e', c.dt);

      c.originalDynamicPower = dlmread(c.powerTrace, '', 1, 0)';

      if ~isempty(c.steps)
        c.dynamicPower = Utils.replicate(c.originalDynamicPower, c.steps);
      else
        c.dynamicPower = c.originalDynamicPower;
        c.steps = size(c.dynamicPower, 2);
      end

      assert(size(c.dynamicPower, 1) == c.cores, ...
        'The number of cores is invalid.');
    end

    function adjustPowerSteps(c, steps)
      c.dynamicPower = Utils.replicate(c.originalDynamicPower, steps);
      c.steps = steps;
    end

    function I = increaseResolution(c, divide)
      dynamicPower = zeros(c.cores, c.steps * divide);

      I = ((1:c.steps) - 1) * divide + 1;
      for i = 1:divide
        dynamicPower(:, I + i - 1) = c.dynamicPower;
      end
      I = I + (divide - 1);

      c.steps = c.steps * divide;
      c.dynamicPower = dynamicPower;

      c.adjustSamplingInterval(c.dt / divide);
    end

    function adjustSamplingInterval(c, dt)
      c.hotspotLine = sprintf('sampling_intvl %.2e', dt);
      c.dt = dt;
    end

    function display(c)
      fprintf('Number of cores:          %d\n', c.cores);
      fprintf('Number of steps:          %d\n', c.steps);
      fprintf('Sampling interval:        %.2e\n', c.dt);
      fprintf('Total simulated time:     %.2f s\n', c.dt * c.steps);
      fprintf('Number of samples for MC: %d\n', c.samples);
      fprintf('Polynomial order for PC:  %d\n', c.order);
    end

    function s = stamp(c, varargin)
      p = struct(...
        'prefix', [], ...
        'cores', c.cores, ...
        'steps', c.steps, ...
        'dt', c.dt);

      for i = 1:floor(length(varargin) / 2)
        name = lower(varargin{(i - 1) * 2 + 1});
        value = varargin{(i - 1) * 2 + 2};
        p.(name) = value;
      end

      s = sprintf('nc%d_ns%d_f%d', p.cores, p.steps, 1 / p.dt);
      if ~isempty(p.prefix), s = [ p.prefix, '_', s ]; end
    end

    function time = timeLine(c)
      time = (0:(c.steps - 1)) * c.dt;
    end

    function set = hotspotSet(c)
      set = { c.floorplan, c.hotspotConfig, c.hotspotLine };
    end
  end
end
