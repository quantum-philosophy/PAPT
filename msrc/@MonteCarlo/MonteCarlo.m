classdef MonteCarlo < handle
  methods (Static)
    [ E, V, out ] = sample(f, dims, samples);
    [ E, V, out, t ] = sample3D(f, dims, samples, stamp);
  end
end
