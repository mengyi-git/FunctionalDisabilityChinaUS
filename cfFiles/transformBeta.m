function y = transformBeta(coef)

% coef in Hanewald et al. (2019) see Table 3
y = coef;

% 3 5 6 
y(3, :) = coef(3, :)/1e2;
y(5, :) = coef(5, :)/1e2;
y(6, :) = coef(6, :)/1e5;

y = y';