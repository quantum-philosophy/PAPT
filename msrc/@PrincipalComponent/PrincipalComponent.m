classdef PrincipalComponent < handle
  methods (Static)
    function [ sdim, map ] = perform(ddim)
      sdim = ddim + 1;

      if ddim == sdim
        %
        % All r.v.'s are independent.
        %
        map = diag(ones(ddim, 1));
      elseif ddim + 1 == sdim
        %
        % The last r.v. is shared, the rest are independent.
        % The ratio is 1:1.
        %
        map = [ diag(ones(ddim, 1) * 0.5) (ones(ddim, 1) * 0.5) ];
      else
        error('Not supported.');
      end
    end
  end
end
