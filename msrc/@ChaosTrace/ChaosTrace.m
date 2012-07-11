classdef ChaosTrace < handle
  properties
    coeff
  end

  methods
    function ct = ChaosTrace(cores, terms, steps)
      ct.coeff = zeros(cores, terms, steps);
    end
  end
end
