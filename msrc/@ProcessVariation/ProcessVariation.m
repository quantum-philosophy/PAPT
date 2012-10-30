classdef ProcessVariation < handle
  properties (Constant)
    %
    % The standard deviation of variations.
    %
    deviationRatio = 0.05;

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
    % variations. (It is good if they sum up to one.)
    %
    localRatio = 0.5;
    globalRatio = 0.5;
  end

  methods (Static)
    function [ P, dimension ] = analyze(floorplan, nominal, varargin)
      %
      % Description:
      %
      %   Does the same as `performPCA', but on the correlation
      %   matrix of the given floorplan.
      %
      C = ProcessVariation.computeCorrelation(floorplan);

      lCount = size(C, 1);
      gCount = ProcessVariation.globalCount;

      deviation = ProcessVariation.deviationRatio * nominal;

      %
      % Take into account the proportions between the two parts.
      %
      lDeviation = sqrt(ProcessVariation.localRatio) * deviation;
      gDeviation = sqrt(ProcessVariation.globalRatio) * deviation;

      S = ones(1, lCount) * lDeviation;

      P = ProcessVariation.performPCA(diag(S) * C * diag(S), varargin{:});

      %
      % Append the global r.v.'s.
      %
      P((end + 1):(end + gCount), (end + 1):(end + gCount)) = ...
        diag(ones(1, gCount) * gDeviation);

      %
      % The last part is mapping of the r.v.'s to the cores.
      %
      M = [ diag(ones(1, lCount)), ones(lCount, gCount) ];
      P = M * P;

      dimension = size(P, 2);
    end
  end

  methods (Static, Access = 'private')
    function P = performPCA(M, reduction)
      %
      % Description:
      %
      %   Performs the PCA on the given matrix and shrinks
      %   the number of principal components according to
      %   the threshold, which is constant for now.
      %
      [ P, L, E ] = pcacov(M);

      if nargin < 2 || isempty(reduction), reduction = 'adjustable'; end

      switch lower(reduction)
      case 'adjustable'
        keep = min(find((cumsum(E) - ProcessVariation.threshold) > 0));
        if isempty(keep), keep = size(P, 1); end
      case 'none'
        keep = size(P, 1);
      otherwise
        error('The reduction method is unknown.');
      end

      P = P(:, 1:keep) * diag(sqrt(L(1:keep)));
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
          C(i, j) = exp(-abs(distance(i) - distance(j)) / corrLength);
          C(j, i) = C(i, j);
        end
      end
    end
  end
end
