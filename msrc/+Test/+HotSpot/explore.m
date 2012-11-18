function explore
  clear all;
  setup;

  use('SystemSimulation');

  processorCount = 2;
  taskCount = 40;

  options = configure('processorCount', processorCount);

  [ platform, application ] = parseTGFF( ...
    File.join(File.trace, '..', 'Assets', ...
      sprintf('%03d_%03d.tgff', processorCount, taskCount)));

  schedule = Schedule.Dense(platform, application);

  hotspot = HotSpot.PowerStepwiseChaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, options.chaosOptions);

  display(hotspot);

  power = DynamicPower(options.samplingInterval);

  function drawSchedule(schedule, title)
    Pdyn = options.powerScale * power.compute(schedule);
    time = options.samplingInterval * (0:(size(Pdyn, 2) - 1));
    [ Texp, Tvar, ~, Pexp, ~, ~ ] = hotspot.computeWithLeakage(Pdyn, options.leakage);

    Utils.drawTemperature(time, Utils.toCelsius(Texp), Tvar, 'layout', 'joint');
    Plot.title('%s: Temperature profile', title);

    power.plot(Pexp);
    Plot.title('%s: Power profile', title);
    Plot.label('Time, s', 'Power, W');
  end

  drawSchedule(schedule, 'Initial');

  populationSize = 10;
  mutationRate = 0.01;

  function population = populate(Genomelength, FitnessFcn, options_)
    M = randi(processorCount, populationSize, taskCount);
    P = rand(populationSize, taskCount);
    population = [ M, P ];
  end

  function children = mutate(parents, options_, nvars, ...
    FitnessFcn, state, thisScore, thisPopulation)

    count = length(parents);
    children = zeros(count, 2 * taskCount);
    for i = 1:count
      child = thisPopulation(parents(i), :);

      %
      % Mutate the mapping part.
      %
      points = find(rand(1, taskCount) < mutationRate);
      child(points) = randi(processorCount, 1, length(points));
      %
      % Mutate the priority part.
      %
      points = find(rand(1, taskCount) < mutationRate);
      child(taskCount + points) = rand(1, length(points));

      children(i, :) = child;
    end
  end

  function fitness = evaluate(chromosome)
    schedule = Schedule.Dense(platform, application, ...
      'mapping', chromosome(1:taskCount), ...
      'priority', chromosome((taskCount + 1):end));
    Pdyn = options.powerScale * power.compute(schedule);
    [ Texp, ~, ~, Pexp, ~, ~ ] = hotspot.computeWithLeakage(Pdyn, options.leakage);
    fitness = sum(Pexp(:)) * options.samplingInterval;
  end

  gaOptions = gaoptimset;
  gaOptions.CreationFcn = @populate;
  gaOptions.CrossoverFcn = @crossoversinglepoint;
  gaOptions.CrossoverFraction = 0.8;
  gaOptions.Display = 'diagnose';
  gaOptions.EliteCount = floor(0.05 * populationSize);
  gaOptions.Generations = 10;
  gaOptions.MigrationFraction = 0.2;
  gaOptions.MutationFcn = @mutate;
  gaOptions.PopulationSize = populationSize;
  gaOptions.SelectionFcn = @selectiontournament;

  tic;
  best = ga(@evaluate, 2 * taskCount, gaOptions);
  fprintf('Genetic algorithm: %.2f s\n', toc);

  schedule = Schedule.Dense(platform, application, ...
    'mapping', best(1:taskCount), ...
    'priority', best((taskCount + 1):end));

  drawSchedule(schedule, 'Optimized');
end
