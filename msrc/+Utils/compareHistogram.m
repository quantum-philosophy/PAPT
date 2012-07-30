function compareHistogram(outMC, outPC, labels)
  [ samplesMC, ddim ] = size(outMC);
  [ samplesPC, ~ ] = size(outPC);

  figure;

  for i = 1:ddim
    p = subplot(1, ddim, i);
    title('Histogram');

    mc = outMC(:, i);
    pc = outPC(:, i);

    x = Utils.constructLinearSpace(mc, pc);

    mcHist = histc(mc, x) / samplesMC;
    pcHist = histc(pc, x) / samplesPC;

    c = Utils.pickColor(1);
    bar(x, mcHist, 'FaceColor', c, 'Edgecolor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    hold on;

    c = Utils.pickColor(2);
    bar(x, pcHist, 'FaceColor', c, 'EdgeColor', c);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.75);

    legend(labels{:});
  end
end
