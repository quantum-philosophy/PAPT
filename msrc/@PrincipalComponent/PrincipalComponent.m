classdef PrincipalComponent < handle
  properties (Constant)
    %
    % This constant determines the percentage of information that
    % kept principal components should explain.
    %
    threshold = 99;

    %
    % The number of global r.v.'s.
    %
    globalCount = 1;

    %
    % The ratio of the global and local variations in the total
    % variations. (It is good when they sum up to 1.)
    %
    globalRatio = 0.5;
    localRatio = 0.5;
  end

  methods (Static)
    function [ P, sdim ] = perform(floorplan, Ldev)
      %
      % Description:
      %
      %   Does the same as `perform', but on the correlation
      %   matrix of the given floorplan.
      %
      C = PrincipalComponent.computeCorrelation(floorplan);

      lCount = size(C, 1);
      gCount = PrincipalComponent.globalCount;

      %
      % Take into account the global variations.
      %
      C((end + 1):(end + gCount), (end + 1):(end + gCount)) = diag(ones(1, gCount));

      if nargin < 2, Ldev = Leakage.Base.Ldev; end

      %
      % Take into account the proportions between the two parts.
      %
      lLdev = PrincipalComponent.localRatio * Ldev;
      gLdev = PrincipalComponent.globalRatio * Ldev;

      S = [ ones(1, lCount) * lLdev, ones(1, gCount) * gLdev ];

      P = PrincipalComponent.performRegular(diag(S) * C * diag(S));

      %
      % The last part is mapping the r.v.'s to the cores.
      %
      M = [ diag(ones(1, lCount)), ones(lCount, gCount) ];
      P = M * P;

      sdim = size(P, 2);
    end
  end

  methods (Static, Access = 'private')
    function P = performRegular(M)
      %
      % Description:
      %
      %   Performs the PCA on the given matrix and shrinks
      %   the number of principal components according to
      %   the threshold, which is constant for now.
      %
      [ P, L, E ] = pcacov(M);
      toKeep = min(find((cumsum(E) - PrincipalComponent.threshold) > 0));
      P = P(:, 1:toKeep) * diag(sqrt(L(1:toKeep)));
    end

    function C = computeCorrelation(floorplan)
      D = dlmread(floorplan, '', 0, 1);

      W = D(:, 1);
      H = D(:, 2);
      X = D(:, 3);
      Y = D(:, 4);

      cores = size(D, 1);

      dieW = max(X + W);
      dieH = max(Y + H);

      dieX = dieW / 2;
      dieY = dieH / 2;

      coreX = X + W / 2;
      coreY = Y + H / 2;

      corrLength = max(dieW, dieH) / 2;

      distance = sqrt((dieX - coreX).^2 + (dieY - coreY).^2);

      C = diag(ones(1, cores));

      for i = 1:cores
        for j = (i + 1):cores
          C(i, j) = exp(- abs(distance(i) - distance(j)) / corrLength);
          C(j, i) = C(i, j);
        end
      end
    end
  end
end
