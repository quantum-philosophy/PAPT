function [ globalError, localError ] = comparePDF(outMC, outPC, time, labels, varargin)
  [ samplesMC, ddim, tdim ] = size(outMC);
  [ samplesPC, ddimPC, tdimPC ] = size(outPC);

  assert(ddim == ddimPC, 'The deterministic dimension is invalid.');
  assert(tdim == tdimPC, 'The time dimension is invalid.');

  draw = false;

  if nargin > 2 && ~isempty(time)
    draw = true;

    if nargin < 4
      labels = {};
      for i = 1:ddim
        labels{end + 1} = sprintf('PE %d', i);
      end
    end

    figure;
    title('Probability Density Error');
  end

  done = 0;
  total = ddim * tdim;

  h = ibar('Comparison of PDFs: step %d out of %d.', total);

  localError = zeros(ddim, tdim);

  for i = 1:ddim
    for j = 1:tdim
      mc = outMC(:, i, j);
      pc = outPC(:, i, j);

      x = Utils.constructLinearSpace(mc, pc, varargin{:});

      densityMC = histc(mc, x) / samplesMC;
      densityPC = histc(pc, x) / samplesPC;

      % densityMC = ksdensity(mc, x);
      % densityPC = ksdensity(pc, x);

      localError(i, j) = Utils.NRMSE(densityMC, densityPC);

      done = done + 1;
      increase(h);
    end

    if draw
      color = Utils.pickColor(i);
      line(time, localError(i, :) * 100, 'Color', color);
    end
  end

  globalError = sqrt(sum(localError(:) .^ 2) / (ddim * tdim));

  if draw
    legend(labels{:});
    xlabel('Time, s');
    ylabel('NRMSE, %');
  end
end
