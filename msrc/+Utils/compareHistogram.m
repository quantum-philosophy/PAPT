function compareHistogram(mcRaw, pcRaw, labels)
  [ mcCount, ddim ] = size(mcRaw);
  [ pcCount, ~ ] = size(pcRaw);

  figure;

  for i = 1:ddim
    p = subplot(1, ddim, i);

    mcraw = mcRaw(:, i);
    pcraw = pcRaw(:, i);

    [ x, dx ] = Utils.constructLinearSpace(mcraw, pcraw);

    mcHist = histc(mcraw, x);
    pcHist = histc(pcraw, x);

    mcDensity = mcHist / (mcCount * dx);
    pcDensity = pcHist / (pcCount * dx);

    c = Utils.pickColor(1);
    bar(x, mcDensity, 'FaceColor', c, 'Edgecolor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    hold on;

    c = Utils.pickColor(2);
    bar(x, pcDensity, 'FaceColor', c, 'EdgeColor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    error = Utils.NRMSE(mcDensity, pcDensity);

    title('Empirical PDF');
    xlabel('Temperature, C');
    ylabel('Empirical PDF');

    if nargin < 3, continue; end

    labelPC = sprintf('%s (NRMSE %.2e)', labels{2}, error);
    legend(labels{1}, labelPC);
  end
end
