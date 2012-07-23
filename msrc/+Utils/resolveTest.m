function [ floorplan, powerTrace, config, configLine ] = resolveTest(cores)
  stamp = sprintf('%02d', cores);

  floorplan = Utils.resolvePath([ stamp, '.flp' ], 'test');
  powerTrace = Utils.resolvePath([ stamp, '.ptrace' ], 'test');
  config = Utils.resolvePath('hotspot.config');
  configLine = 'sampling_intvl 1e-4';
end
