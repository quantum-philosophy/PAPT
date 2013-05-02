function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, 'Assets');

  %
  % Application
  %
  processorCount = options.getSet('processorCount', 4);
  taskCount = options.getSet('taskCount', 20 * processorCount);

  tgffConfig = File.join('+Test', 'Assets', ...
    sprintf('%03d_%03d.tgff', processorCount, taskCount));

  [ options.platform, options.application ] = Utils.parseTGFF(tgffConfig);
  options.schedule = Schedule.Dense(options.platform, options.application);

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
  % Temperature
  %
  options.temperatureOptions = Options( ...
    'floorplan', File.join(path, sprintf('%03d.flp', processorCount)), ...
    'config', File.join(path, 'hotspot.config'), ...
    'line', sprintf('sampling_intvl %.4e', options.samplingInterval), ...
    'leakage', LeakagePower.LinearInterpolation( ...
      'dynamicPower', options.dynamicPower, ...
      'filename', File.join(path, 'inverter_45nm_L5_T1000.leak')));

  %
  % Polynomil chaos
  %
  options.chaosOptions = Options('order', 4, ...
    'quadratureOptions', Options( ...
      'method', 'adaptive', ...
      'ruleName', 'GaussHermiteHW', ...
      'polynomialOrder', 4));
end
