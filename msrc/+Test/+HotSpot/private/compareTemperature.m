function draw(time, expectationSet, varianceSet)
  setCount = length(expectationSet);
  processorCount = size(expectationSet{1}, 1);

  for i = 1:processorCount
    figure;
    for j = 1:setCount
      color = Color.pick(j);
      line(time, expectationSet{j}(i, :), ...
        'Color', color, 'LineWidth', 1);
      line(time, expectationSet{j}(i, :) + sqrt(varianceSet{j}(i, :)), ...
        'Color', color, 'LineStyle', '--');
    end
    Plot.title('Temperature (PE%d)', i);
    Plot.label('Time, s', 'Temperature, C');
    Plot.limit(time);
  end
end

