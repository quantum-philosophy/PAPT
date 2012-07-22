classdef Constants < handle
  properties (Constant)
    workingDirectory = [ Constants.thisDirectory, '/../build' ];
    vendorDirectory = [ Constants.thisDirectory, '/../vendor' ];
    cacheDirectory = [ Constants.thisDirectory, '/../build/cache' ];
    testDirectory = [ Constants.thisDirectory, '/../test' ];

    zeroKelvin = -273.15; % in degree Celsium

    %
    % Visualization
    %
    roundRobinColors = { ...
      [ 87, 181, 232] / 255, ...
      [ 230, 158, 0 ] / 255, ...
      [ 129, 197, 122 ] / 255, ...
      [ 20, 43, 140 ] / 255, ...
      [ 195, 0, 191 ] / 255, ...
      'c', 'y' };
  end

  methods (Static)
    function result = thisDirectory
      filename = mfilename('fullpath');
      attrs = regexp(filename, '^(.*)/[^/]+$', 'tokens');
      result = attrs{1}{1};
    end
  end
end
