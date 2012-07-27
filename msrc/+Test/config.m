classdef config < handle
  properties
    cores
    steps
    samples
    dt
    floorplan
    powerTrace
    hotspotConfig
    hotspotLine
    dynamicPower
  end

  methods
    function c = config()
      c.cores = 4;
      c.steps = 100;
      c.samples = 1e4;
      c.dt = 1e-3;

      [ c.floorplan, c.powerTrace, c.hotspotConfig ] = ...
        Utils.resolveTest(c.cores);

      c.hotspotLine = sprintf('sampling_intvl %.2e', c.dt);

      c.dynamicPower = dlmread(c.powerTrace, '', 1, 0)';
      c.steps = size(c.dynamicPower, 2);

      assert(size(c.dynamicPower, 1) == c.cores, ...
        'The number of cores is invalid.');
    end

    function adjustPowerSteps(c, steps)
      c.dynamicPower = Utils.replicate(c.dynamicPower, steps);
      c.steps = steps;
    end

    function adjustSamplingInterval(c, dt)
      c.hotspotLine = sprintf('sampling_intvl %.2e', dt);
      c.dt = dt;
    end

    function display(c)
      fprintf('Number of cores:   %d\n', c.cores);
      fprintf('Number of steps:   %d\n', c.steps);
      fprintf('Number of samples: %d\n', c.samples);
      fprintf('Sampling interval: %.2e\n', c.dt);
    end

    function s = stamp(c, prefix)
      if nargin < 2
        s = sprintf('nc%d_ns%d_f%d', c.cores, c.steps, 1 / c.dt);
      else
        s = sprintf('%s_nc%d_ns%d_f%d', prefix, c.cores, c.steps, 1 / c.dt);
      end
    end

    function time = timeLine(c)
      time = (0:(c.steps - 1)) * c.dt;
    end
  end
end
