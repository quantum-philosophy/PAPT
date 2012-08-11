function value = extract(options, field, default)
  if isfield(options, field)
    value = options.(field);
  else
    value = default;
  end
end
