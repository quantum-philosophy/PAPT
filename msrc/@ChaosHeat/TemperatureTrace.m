classdef TemperatureTrace < handle
  properties (SetAccess = 'private')
    %
    % All the coefficients for all the cores for all the steps.
    %
    coeff
  end

  methods
    function tt = TemperatureTrace(steps, cores, terms)
      coeff = zeros(steps, cores, terms);
    end
  end
end
