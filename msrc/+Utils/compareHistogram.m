function compareHistogram(mcRaw, pcRaw, labels)
  [ samplesMC, ddim ] = size(mcRaw);
  [ samplesPC, ~ ] = size(pcRaw);

  figure;

  for i = 1:ddim
    p = subplot(1, ddim, i);
    title('Empirical PDF');

    mcraw = mcRaw(:, i);
    pcraw = pcRaw(:, i);

    x = Utils.constructLinearSpace(mcraw, pcraw);

    mcHist = histc(mcraw, x) / samplesMC;
    pcHist = histc(pcraw, x) / samplesPC;

    c = Utils.pickColor(1);
    bar(x, mcHist, 'FaceColor', c, 'Edgecolor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    hold on;

    c = Utils.pickColor(2);
    bar(x, pcHist, 'FaceColor', c, 'EdgeColor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    error = Utils.NRMSE(mcHist, pcHist);

    labelPC = sprintf('%s (NRMSE %.2e)', labels{2}, error);
    legend(labels{1}, labelPC);
    xlabel('Temperature, C');
    ylabel('Empirical PDF');
  end
end
