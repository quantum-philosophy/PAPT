classdef PC < handle
  properties (SetAccess = 'protected')
    dimension
    order
    count
    psi
    norm
  end

  methods
    function pc = PC(dimension, order)
      if nargin < 2, order = 2; end;

      pc.dimension = dimension;
      pc.order = order;

      [ pc.psi, pc.norm, pc.count ] = PC.calculateExpansion(dimension, order);
    end
  end

  methods (Static)
    function count = calculateCount(dimension, order)
      count = factorial(dimension + order) / ...
        (factorial(dimension) * factorial(order));
    end
  end

  methods (Static, Access = 'private')
    psi = construct1D(x, order);
    psi = constructXD(x, count);

    function [ psi, norm, count ] = calculateExpansion(dimension, order)
      filename = [ 'PC_d', num2str(dimension), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename);

      if exist(filename, 'file')
        load(filename);
      else
        count = PC.calculateCount(dimension, order);

        for i = 1:dimension
          x(i) = sym([ 'x', num2str(i) ], 'real');
        end

        psi = PC.constructXD(x, count);

        norm = zeros(1, count);
        norm(1) = 1;
        for i = 2:count
          norm(i) = GQ.integrate(@(s) subs(psi(i) * psi(i), x, s), dimension);
        end

        save(filename, 'psi', 'norm', 'count');
      end
    end
  end
end
