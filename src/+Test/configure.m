function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, 'Assets');

  %
  % Platform and application
  %
  processorCount = options.getSet('processorCount', 4);
  taskCount = options.getSet('taskCount', 20 * processorCount);

  tgffConfig = File.join('+Test', 'Assets', ...
    sprintf('%03d_%03d.tgff', processorCount, taskCount));

  [ options.platform, options.application ] = Utils.parseTGFF(tgffConfig);
  options.schedule = Schedule.Dense(options.platform, options.application);

  options.die = Die('floorplan', ...
    File.join(path, sprintf('%03d.flp', processorCount)));

  %
  % Dynamic power
  %
  options.power = DynamicPower(options.getSet('samplingInterval', 1e-3));
  options.dynamicPower = options.getSet('powerScale', 1) * ...
    options.power.compute(options.schedule);

  if options.has('stepCount')
    options.dynamicPower = Utils.stretch( ...
      options.dynamicPower, options.stepCount);
  else
    options.stepCount = size(options.dynamicPower, 2);
  end

  resample = options.get('resample', 1);
  if resample > 1
    options.stepCount = options.stepCount * resample;
    options.samplingInterval = options.samplingInterval / resample;
    options.dynamicPower = Utils.resample(options.dynamicPower, resample);
  end

  %
  % Leakage power
  %
  options.leakage = LeakagePower.LinearInterpolation( ...
    'dynamicPower', options.dynamicPower, ...
    'filename', File.join(path, 'inverter_45nm_L5_T1000_08.leak'));

  %
  % Process variation
  %
  eta = 0.00;
  lse = 1.00 * options.die.radius;
  lou = 1.00 * options.die.radius;

  function K = correlate(s, t)
    %
    % Squared exponential kernel
    %
    Kse = exp(-sum((s - t).^2, 1) / lse^2);

    %
    % Ornstein-Uhlenbeck kernel
    %
    rs = sqrt(sum(s.^2, 1));
    rt = sqrt(sum(t.^2, 1));
    Kou = exp(-abs(rs - rt) / lou);

    K = eta * Kse + (1 - eta) * Kou;
  end

  options.processModel = 'Normal';
  options.processOptions = Options( ...
    'die', options.die, ...
    'expectation', options.leakage.Lnom, ...
    'deviation', 0.05 * options.leakage.Lnom, ...
    'kernel', @correlate, ...
    'globalPortion', 0.5, ...
    'threshold', 0.99);

  %
  % Temperature
  %
  options.temperatureOptions = Options('die', options.die, ...
    'config', File.join(path, 'hotspot.config'), ...
    'line', sprintf('sampling_intvl %.4e', options.samplingInterval), ...
    'leakage', options.leakage);

  %
  % Polynomil chaos
  %
  options.quadratureOptions = Options( ...
    'method', 'adaptive', ...
    'ruleName', 'GaussHermiteHW', ...
    'polynomialOrder', 4);

  options.chaosOptions = Options('order', 4, ...
    'quadratureOptions', options.quadratureOptions);
end
