function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, 'Assets');

  processorCount = options.getSet('processorCount', 4);

  powerScale = options.getSet('powerScale', 2);

  powerProfile = File.join(path, sprintf('%02d.ptrace', processorCount));
  options.powerProfile = powerScale * dlmread(powerProfile, '', 1, 0).';
  options.processorCount = size(options.powerProfile, 1);

  leakage = File.join(path, 'inverter_45nm.leak');
  options.leakage = LeakagePower( ...
    'dynamicPower', options.powerProfile, ...
    'filename', leakage, ...
    'order', [ 1, 2 ], ...
    'scale', [ 1, 0.7, 0; 1, 1, 1 ]);

  if options.has('stepCount')
    options.powerProfile = Utils.stretch( ...
      options.powerProfile, options.stepCount);
  else
    options.set('stepCount', size(options.powerProfile, 2));
  end

  options.samplingInterval = 1e-3;

  resample = options.get('resample', 1);
  if resample > 1
    options.stepCount = options.stepCount * resample;
    options.samplingInterval = options.samplingInterval / resample;
    options.powerProfile = Utils.resample(options.powerProfile, resample);
  end

  options.hotspotOptions = Options( ...
    'floorplan', File.join(path, sprintf('%02d.flp', processorCount)), ...
    'config', File.join(path, 'hotspot.config'), ...
    'line', sprintf('sampling_intvl %.4e', options.samplingInterval));

  options.chaosOptions = Options('order', 4, ...
    'quadratureOptions', Options( ...
      'method', 'adaptive', ...
      'ruleName', 'GaussHermiteHW', ...
      'polynomialOrder', 4));
end
