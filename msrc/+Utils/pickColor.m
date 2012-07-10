function color = pickColor(i)
  colors = Constants.roundRobinColors;
  color = colors{mod(i - 1, length(colors)) + 1};
end
