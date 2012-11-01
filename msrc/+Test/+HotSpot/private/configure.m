function options = configure(varargin)
  options = Options(varargin{:});

  path = File.join(File.trace, '..', '..', 'Assets');

  processorCount = options.getSet('processorCount', 4);

  options.set('floorplan', ...
    File.join(path, sprintf('%02d.flp', processorCount)));
  options.set('hotspotConfig', File.join(path, 'hotspot.config'));
  options.set('hotspotLine', 'sampling_intvl 1e-3');

  powerProfile = File.join(path, sprintf('%02d.ptrace', processorCount));
  options.set('powerProfile', 2 * dlmread(powerProfile, '', 1, 0).');
  options.set('processorCount', size(options.powerProfile, 1));

  if options.has('stepCount')
    options.powerProfile = Utils.stretch( ...
      options.powerProfile, options.stepCount);
  else
    options.set('stepCount', size(options.powerProfile, 2));
  end

  leakage = File.join(path, 'inverter_45nm.leak');
  options.set('leakage', LeakagePower(options.powerProfile, ...
    'filename', leakage, ...
    'order', [ 1, 2 ], ...
    'scale', [ 1, 0.7, 0; 1, 1, 1 ]));

  options.set('chaosOptions', Options('order', 4, ...
    'quadratureOptions', Options( ...
      'method', 'adaptive', ...
      'ruleName', 'GaussHermiteHW', ...
      'polynomialOrder', 4)));
end
