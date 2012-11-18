function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, '..', '..', 'Assets');

  processorCount = options.getSet('processorCount', 4);

  powerScale = options.getSet('powerScale', 2);

  powerProfile = File.join(path, sprintf('%02d.ptrace', processorCount));
  options.set('powerProfile', powerScale * dlmread(powerProfile, '', 1, 0).');
  options.set('processorCount', size(options.powerProfile, 1));

  leakage = File.join(path, 'inverter_45nm.leak');
  options.set('leakage', LeakagePower(options.powerProfile, ...
    'filename', leakage, ...
    'order', [ 1, 2 ], ...
    'scale', [ 1, 0.7, 0; 1, 1, 1 ]));

  if options.has('stepCount')
    options.powerProfile = Utils.stretch( ...
      options.powerProfile, options.stepCount);
  else
    options.set('stepCount', size(options.powerProfile, 2));
  end

  options.set('samplingInterval', 1e-3);

  resample = options.get('resample', 1);
  if resample > 1
    options.set('stepCount', options.stepCount * resample);
    options.set('samplingInterval', options.samplingInterval / resample);
    options.set('powerProfile', Utils.resample(options.powerProfile, resample));
  end

  options.set('floorplan', ...
    File.join(path, sprintf('%02d.flp', processorCount)));
  options.set('hotspotConfig', File.join(path, 'hotspot.config'));
  options.set('hotspotLine', ...
    sprintf('sampling_intvl %.4e', options.samplingInterval));

  options.set('chaosOptions', Options('order', 5, ...
    'quadratureOptions', Options( ...
      'method', 'adaptive', ...
      'ruleName', 'GaussHermiteHW', ...
      'polynomialOrder', 5)));
end
