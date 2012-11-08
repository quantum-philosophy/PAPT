function scaleProcessors
  compareSpeed( ...
    'name', 'Number of processing elements', ...
    'shortName', 'Processors', ...
    'range', [ 2 4 8 16 32 ], ...
    'repeat', [ 4 4 1 1 1 ], ...
    'configure', @(processorCount) configure( ...
      'processorCount', processorCount, 'stepCount', 1e3));
end