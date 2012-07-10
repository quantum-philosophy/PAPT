classdef Trace < handle
  properties (SetAccess = 'private')
    %
    % Steps, cores, terms, just in case.
    %
    steps
    cores
    terms

    %
    % The leakage power model.
    %
    leakage

    %
    % All the coefficients for all the cores for all the steps.
    %
    coeff

    %
    % The current position in discrete time.
    %
    step
  end

  methods
    function trace = Trace(steps, cores, terms, leakage)
      trace.steps = steps;
      trace.cores = cores;
      trace.terms = terms;
      trace.leakage = leakage;

      trace.coeff = zeros(steps, cores, terms);

      trace.step = 0;
    end

    function advance(trace)
      trace.step = trace.step + 1;
    end

    function P = computeLeakagePower(trace, values)

    end
  end
end
