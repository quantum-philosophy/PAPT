function s = methodStamp(method)
  s = sprintf('PolynomialChaos(%s_co%d)', method.chaosName, method.chaosOrder);
end
