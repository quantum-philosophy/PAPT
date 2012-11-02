function drawTemperature(time, expectationSet, varianceSet, varargin)
  options = Options(varargin{:});

  if ~isa(expectationSet, 'cell')
    expectationSet = { expectationSet };
    varianceSet = { varianceSet };
  end

  switch options.get('method', 'separate')
  case 'separate'
    drawSeparate(time, expectationSet, varianceSet);
  case 'joint'
    drawJoint(time, expectationSet, varianceSet);
  end
end

function drawSeparate(time, expectationSet, varianceSet)
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

function drawJoint(time, expectationSet, varianceSet)
  setCount = length(expectationSet);
  processorCount = size(expectationSet{1}, 1);

  figure;
  for i = 1:processorCount
    color = Color.pick(i);
    for j = 1:setCount
      line(time, expectationSet{j}(i, :), ...
        'Color', color, 'LineWidth', 1);
      line(time, expectationSet{j}(i, :) + sqrt(varianceSet{j}(i, :)), ...
        'Color', color, 'LineStyle', '--');
    end
  end
  Plot.title('Temperature', i);
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(time);
end
