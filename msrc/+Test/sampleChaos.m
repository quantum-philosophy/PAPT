function [ expAn, varAn, raw, expEm, varEm ] = sampleChaos(ch, c)
  [ expAn, varAn, trace ] = ch.solve(c.dynamicPower);
  expAn = Utils.toCelsius(expAn);

  if nargout < 3, return; end

  raw = ch.pc.evaluateSet(trace, normrnd(0, 1, ch.pc.sdim, c.polynomialChaosSamples));
  raw = Utils.toCelsius(raw);

  if nargout < 4, return; end

  expEm = zeros(size(expAn));
  varEm = zeros(size(varAn));

  for i = 1:c.steps
    snapshot = raw(:, :, i);
    expEm(:, i) = mean(snapshot);
    varEm(:, i) = var(snapshot);
  end
end
