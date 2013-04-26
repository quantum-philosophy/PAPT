function [ T, P ] = solveDynamicSteadyStateWithLeakage(this, Pdyn, options)
  nodeCount = this.nodeCount;
  [ processorCount, stepCount ] = size(Pdyn);
  assert(processorCount == this.processorCount);

  E = this.E;
  D = this.D;
  BT = this.BT;
  U = this.U;
  UT = this.UT;
  Lambda = this.L;

  Tamb = this.ambientTemperature;
  dt = this.samplingInterval;

  leak = options.leakage.evaluate;
  L = options.L;

  iterationLimit = options.get('iterationLimit', 10);
  tolerance = options.get('tolerance', 0.5);

  function T_ = computeOne(P_)
    W = zeros(nodeCount, stepCount);
    X = zeros(nodeCount, stepCount);

    Q = D * P_;

    W(:, 1) = Q(:, 1);

    for i = 2:stepCount
      W(:, i) = E * W(:, i - 1) + Q(:, i);
    end

    X(:, 1) = U * diag(1 ./ (1 - exp(dt * ...
      stepCount * Lambda))) * UT * W(:, stepCount);

    for i = 2:stepCount
      X(:, i) = E * X(:, i - 1) + Q(:, i - 1);
    end

    T_ = BT * X + Tamb;
  end

  sampleCount = size(L, 2);

  T = zeros(processorCount * stepCount, sampleCount);
  P = zeros(processorCount * stepCount, sampleCount);

  %
  % NOTE: 'i' is already occupied!
  %
  for j = 1:sampleCount
    l = Utils.replicate(L(:, j), 1, stepCount);

    Pcurrent = Pdyn + leak(l, Tamb);
    Tcurrent = computeOne(Pcurrent);

    for k = 2:iterationLimit
      Tlast = Tcurrent;

      Pcurrent = Pdyn + leak(l, Tcurrent);
      Tcurrent = computeOne(Pcurrent);

      if max(abs(Tcurrent(:) - Tlast(:))) < tolerance, break; end
    end

    T(:, j) = Tcurrent(:);
    P(:, j) = Pcurrent(:);
  end
end
