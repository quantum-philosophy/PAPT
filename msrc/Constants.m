classdef Constants < handle
  properties (Constant)
    workingDirectory = [ Constants.thisDirectory, '/../build' ];

    zeroKelvin = -273.15; % in degree Celsium
  end

  methods (Static)
    function result = thisDirectory
      filename = mfilename('fullpath');
      attrs = regexp(filename, '^(.*)/[^/]+$', 'tokens');
      result = attrs{1}{1};
    end
  end
end
