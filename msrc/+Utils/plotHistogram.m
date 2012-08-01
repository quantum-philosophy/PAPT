function plotHistogram(Raw, varargin)
  [ samples, ddim ] = size(Raw);

  figure;

  color = Utils.pickColor(1);

  for i = 1:ddim
    subplot(1, ddim, i);
    title('Empirical PDF');

    raw = Raw(:, i);

    x = Utils.constructLinearSpace(raw, varargin{:});

    hist = histc(raw, x) / samples;

    bar(x, hist, 'FaceColor', color, 'Edgecolor', 'w');

    xlabel('Temperature, C');
    ylabel('Empirical PDF');
  end
end
