function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, 'Assets');

  %
  % Application
  %
  processorCount = options.getSet('processorCount', 4);

  tgffConfig = File.join('+Test', 'Assets', ...
    sprintf('%03d_%03d.tgff', processorCount, 20 * processorCount));

  [ platform, application ] = parseTGFF(tgffConfig);
  options.schedule = Schedule.Dense(platform, application);

  %
  % Dynamic power
  %
  options.power = DynamicPower(options.getSet('samplingInterval', 1e-3));
  options.powerProfile = options.getSet('powerScale', 1) * ...
    options.power.compute(options.schedule);

  if options.has('stepCount')
    options.powerProfile = Utils.stretch( ...
      options.powerProfile, options.stepCount);
  else
    options.stepCount = size(options.powerProfile, 2);
  end

  resample = options.get('resample', 1);
  if resample > 1
    options.stepCount = options.stepCount * resample;
    options.samplingInterval = options.samplingInterval / resample;
    options.powerProfile = Utils.resample(options.powerProfile, resample);
  end

  %
  % Leakage power
  %
  options.leakage = LeakagePower( ...
    'dynamicPower', options.powerProfile, ...
    'filename', File.join(path, 'inverter_45nm.leak'), ...
    'order', [ 1, 2 ], ...
    'scale', [ 1, 0.7, 0; 1, 1, 1 ]);

  %
  % Temperature
  %
  options.hotspotOptions = Options( ...
    'floorplan', File.join(path, sprintf('%03d.flp', processorCount)), ...
    'config', File.join(path, 'hotspot.config'), ...
    'line', sprintf('sampling_intvl %.4e', options.samplingInterval));

  %
  % Polynomil chaos
  %
  options.chaosOptions = Options('order', 4, ...
    'quadratureOptions', Options( ...
      'method', 'adaptive', ...
      'ruleName', 'GaussHermiteHW', ...
      'polynomialOrder', 4));
end
