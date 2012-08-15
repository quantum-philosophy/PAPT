function observe(varargin)
  [ raw, options ] = Options.extract(varargin{:});
  assert(length(raw) == 1, 'The observation is supported only for one set of data.');
  observe2D(raw{1}, options);
end

function observe2D(Raw, options)
  [ ~, ddim ] = size(Raw);

  figure;

  for i = 1:ddim
    p = subplot(1, ddim, i);

    raw = Raw(:, i);

    x = Stats.constructLinearSpace(raw, options);
    [ x, data ] = Stats.process(x, raw, options);

    Stats.draw(x, data, options);
  end
end
