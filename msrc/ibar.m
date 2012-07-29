classdef ibar < handle
  properties (Access = 'private')
    template
    total
    done
    handle
  end

  methods
    function i = ibar(template, total, done)
      if nargin < 3, done = 0; end
      i.template = template;
      i.total = total;
      i.done = done;

      color = Utils.pickColor(randi(10));

      i.handle = waitbar(done / total, ...
        sprintf(i.template, done, total));
      set(findobj(i.handle, 'type', 'patch'), ...
        'edgecolor', color, 'facecolor', color);
    end

    function increase(i)
      i.done = i.done + 1;
      waitbar(i.done / i.total, i.handle, ...
        sprintf(i.template, i.done, i.total));
      if i.done == i.total, close(i); end
    end

    function close(i)
      close(i.handle);
    end
  end
end
