function [ T, P ] = solveTransientWithLeakage(this, Pdyn, options)
  [ processorCount, stepCount ] = size(Pdyn);
  assert(processorCount == this.processorCount);

  E = this.E;
  D = this.D;
  BT = this.BT;
  Tamb = this.ambientTemperature;

  leak = options.leakage.evaluate;
  L = options.L;

  sampleCount = size(L, 2);

  T = zeros(processorCount * stepCount, sampleCount);
  P = zeros(processorCount * stepCount, sampleCount);

  P_ = bsxfun(@plus, Pdyn(:, 1), leak(L, Tamb));
  X_ = D * P_;
  T_ = BT * X_ + Tamb;

  range = 1:processorCount;
  T(range, :) = T_;
  P(range, :) = P_;

  for i = 2:stepCount
    P_ = bsxfun(@plus, Pdyn(:, i), leak(L, T_));
    X_ = E * X_ + D * P_;
    T_ = BT * X_ + Tamb;

    range = range + processorCount;
    T(range, :) = T_;
    P(range, :) = P_;
  end
end
