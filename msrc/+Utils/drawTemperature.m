function drawTemperature(time, expectationSet, varianceSet, varargin)
  if nargin < 3
    varianceSet = [];
    varargin = {};
  end

  options = Options(varargin{:});

  if ~isa(expectationSet, 'cell')
    expectationSet = { expectationSet };
    varianceSet = { varianceSet };
  end

  switch options.get('layout', 'separate')
  case 'separate'
    drawSeparate(time, expectationSet, varianceSet, options);
  case 'joint'
    drawJoint(time, expectationSet, varianceSet, options);
  end
end

function drawSeparate(time, expectationSet, varianceSet, options)
  setCount = length(expectationSet);
  processorCount = size(expectationSet{1}, 1);

  labels = options.get('labels', cell(1, setCount));

  for i = 1:processorCount
    figure;
    legend = {};
    for j = 1:setCount
      color = Color.pick(j);
      line(time, expectationSet{j}(i, :), ...
        'Color', color, 'LineWidth', 1);

      if isempty(varianceSet{j})
        legend{end + 1} = labels{j};
        continue;
      else
        legend{end + 1} = sprintf('%s: mean', labels{j});
      end

      line(time, expectationSet{j}(i, :) + sqrt(varianceSet{j}(i, :)), ...
        'Color', color, 'LineStyle', '--');
      legend{end + 1} = ...
        sprintf('%s: mean + sigma', labels{j});
    end
    Plot.title('Temperature (PE%d)', i);
    Plot.label('Time, s', 'Temperature, C');
    Plot.limit(time);
    Plot.legend(legend{:});
  end
end

function drawJoint(time, expectationSet, varianceSet, options)
  setCount = length(expectationSet);
  processorCount = size(expectationSet{1}, 1);

  labels = options.get('labels', cell(1, setCount));

  figure;
  legend = {};
  for i = 1:processorCount
    color = Color.pick(i);
    for j = 1:setCount
      line(time, expectationSet{j}(i, :), ...
        'Color', color, 'LineWidth', 1);

      if isempty(varianceSet{j})
        legend{end + 1} = sprintf('%s: PE%d', labels{j}, i);
        continue;
      else
        legend{end + 1} = sprintf('%s: PE%d: mean', labels{j}, i);
      end

      line(time, expectationSet{j}(i, :) + sqrt(varianceSet{j}(i, :)), ...
        'Color', color, 'LineStyle', '--');
      legend{end + 1} = ...
        sprintf('%s: PE%d: mean + sigma (%s)', labels{j}, i);
    end
  end
  Plot.title('Temperature', i);
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(time);
  Plot.legend(legend{:});
end