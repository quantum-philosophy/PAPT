classdef Base < handle
  properties (Constant)
    %
    % The nominal value of the channel length and its deviation.
    %
    Lnom = 45e-9;
    Ldev = 0.05 * Leakage.Base.Lnom;

    %
    % The following constants are used to construct an instance
    % of the leakage model that produces `PleakPdyn' portion
    % of the given dynamic power at temperature `Tref'.
    %
    PleakPdyn = 2/3;
    Tref = Utils.toKelvin(120);
  end

  properties (Access = 'protected')
    %
    % The temperature of the ambience; it is used for the very
    % first sampling round.
    %
    Tamb

    %
    % The mapping matrix from r.v.'s to cores.
    %
    pca

    %
    % One of the leakage parameters.
    %
    alpha
  end

  methods
    function lk = Base(Tamb, cores, pca)
      [ ddim, sdim ] = size(pca);

      assert(ddim == cores, 'The dimensions do not match.');

      lk.Tamb = ones(cores, 1) * Tamb;
      lk.pca = pca;
    end
  end
end
