classdef MonteCarlo < handle
  methods (Static)
    [ E, V, out ] = perform(f, dims, samples);
    [ E, V, t ] = perform3D(f, dims, samples, stamp);
  end
end
