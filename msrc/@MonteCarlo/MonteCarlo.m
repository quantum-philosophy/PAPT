classdef MonteCarlo < handle
  methods (Static)
    [ E, C, out ] = perform(f, dims, samples);
    [ E, C ] = perform3D(f, dims, samples);
  end
end
